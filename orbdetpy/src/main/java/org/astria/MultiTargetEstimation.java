/*
 * MultiTargetEstimation.java - JPDA, CAR-MHF implementation.
 * Copyright (C) 2019-2022 University of Texas
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

package org.astria;

import java.util.Arrays;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import org.hipparchus.distribution.continuous.ChiSquaredDistribution;
import org.hipparchus.geometry.euclidean.threed.Vector3D;
import org.hipparchus.linear.Array2DRowRealMatrix;
import org.hipparchus.linear.ArrayRealVector;
import org.hipparchus.linear.CholeskyDecomposition;
import org.hipparchus.linear.DiagonalMatrix;
import org.hipparchus.linear.MatrixUtils;
import org.hipparchus.linear.RealMatrix;
import org.hipparchus.linear.RealVector;
import org.hipparchus.util.FastMath;
import org.orekit.estimation.measurements.ObservedMeasurement;
import org.orekit.estimation.measurements.Range;
import org.orekit.orbits.CartesianOrbit;
import org.orekit.propagation.SpacecraftState;
import org.orekit.time.AbsoluteDate;
import org.orekit.utils.Constants;
import org.orekit.utils.PVCoordinates;
import org.orekit.utils.TimeStampedPVCoordinates;

public final class MultiTargetEstimation
{
    public static class MultiTargetOutput
    {
        public ArrayList<ArrayList<Estimation.EstimationOutput>> estOutput;
        public ArrayList<ArrayList<Integer>> associatedObs;
        public ArrayList<Integer> unassociatedObs;
    }

    private final Settings odCfg;
    private final ArrayList<Settings> cfgList;
    private final Measurements odObs;

    private final int measSize;
    private final Settings.MeasurementType[] measNames;
    private final boolean singleObject;

    private MultiTargetOutput multiOutput;

    public MultiTargetEstimation(ArrayList<Settings> cfgList, ArrayList<Measurements> obsList)
    {
        this.cfgList = cfgList;
        this.odCfg = cfgList.get(0);
        this.odObs = obsList.get(0);
        measNames = odCfg.cfgMeasurements.keySet().toArray(new Settings.MeasurementType[0]);
        Arrays.sort(measNames);
        switch (measNames[0])
        {
        case AZIMUTH:
        case RANGE:
        case RANGE_RATE:
        case RIGHT_ASCENSION:
            measSize = measNames.length;
            break;
        case POSITION:
            measSize = 3;
            break;
        case POSITION_VELOCITY:
            measSize = 6;
            break;
        default:
            throw(new RuntimeException("Invalid measurement type"));
        }
        singleObject = (measNames[0] != Settings.MeasurementType.RANGE && measNames[0] != Settings.MeasurementType.RANGE_RATE) ||
            measNames.length == 1;
    }

    public MultiTargetOutput multiTargetDetermineOrbit()
    {
        multiOutput = new MultiTargetOutput();
        new MultiTargetFilter().determineOrbit();
        return(multiOutput);
    }

    private class MultiTargetFilter
    {
        private int numStates, numSigmas;
        private HashMap<String, Integer> biasPos;
        private ManualPropagation propagator;

        private class SmootherTimeStep
        {
            RealMatrix Ppre;
            RealMatrix Ppost;
            RealMatrix xpre;
            RealMatrix xpost;
            RealMatrix sigPre;
            RealMatrix sigPost;

            RealMatrix xstar;
            RealMatrix Pstar;

            AbsoluteDate tmSmoother;
            ObservedMeasurement[] measObjs;
        }

        private class Hypothesis
        {
            RealVector xhat;
            RealMatrix P;
            double weight;

            public Hypothesis(RealVector xhat, RealMatrix P, double weight)
            {
                this.xhat = xhat;
                this.P = P;
                this.weight = weight;
            }
        }

        private class JPDALikelihood
        {
            int object;
            int measurement;
            double psi;

            public JPDALikelihood(int object, int measurement, double psi)
            {
                this.object = object;
                this.measurement = measurement;
                this.psi = psi;
            }
        }

        private class SpaceObject
        {
            RealMatrix P;
            Array2DRowRealMatrix sigma;
            Array2DRowRealMatrix propSigma;
            Array2DRowRealMatrix estimMeas;

            RealMatrix Pprop;
            RealVector xhat;
            RealVector xhatPrev;
            double hypothesisWeight;

            ArrayList<SmootherTimeStep> smootherData;
            ArrayList<JPDALikelihood> marginalEvents;
            ArrayList<Integer> associatedObsIndex;
            ArrayList<Estimation.EstimationOutput> estOutput;

            boolean dataAssociated;
            boolean isConsistent;
            boolean fromCAR;

            public SpaceObject(Settings cfg, int Rsize)
            {
                P = new DiagonalMatrix(cfg.estmCovariance);
                sigma = new Array2DRowRealMatrix(numStates, numSigmas);
                propSigma = new Array2DRowRealMatrix(numStates, numSigmas);
                estimMeas = new Array2DRowRealMatrix(Rsize, numSigmas);
                Pprop = null;
                double[] x0 = cfg.getInitialState();
                xhat = new ArrayRealVector(x0);
                xhatPrev = new ArrayRealVector(x0);
                smootherData = new ArrayList<SmootherTimeStep>();
                marginalEvents = new ArrayList<JPDALikelihood>();
                associatedObsIndex = new ArrayList<Integer>();
                estOutput = new ArrayList<Estimation.EstimationOutput>();
                dataAssociated = true;
                hypothesisWeight = 1;
                fromCAR = false;
            }

            public SpaceObject(double[] x, RealMatrix cov, int Rsize, double hypWeight)
            {
                P = cov;
                sigma = new Array2DRowRealMatrix(numStates, numSigmas);
                propSigma = new Array2DRowRealMatrix(numStates, numSigmas);
                estimMeas = new Array2DRowRealMatrix(Rsize, numSigmas);
                Pprop = null;
                xhat = new ArrayRealVector(x);
                xhatPrev = new ArrayRealVector(x);
                smootherData = new ArrayList<SmootherTimeStep>();
                marginalEvents = new ArrayList<JPDALikelihood>();
                associatedObsIndex = new ArrayList<Integer>();
                estOutput = new ArrayList<Estimation.EstimationOutput>();
                dataAssociated = true;
                hypothesisWeight = hypWeight;
                fromCAR = true;
            }
        }

        private void determineOrbit()
        {
            numStates = odCfg.parameters.size() + 6;
            numSigmas = 2*numStates;
            boolean activateCAR = false;
            AbsoluteDate tm = odCfg.propStart;
            final double weight = 0.5/numStates;
            propagator = new ManualPropagation(odCfg);
            final SpacecraftState[] ssta = new SpacecraftState[1];
            final ArrayList<SpaceObject> promotedTracks = new ArrayList<SpaceObject>();
            final ArrayList<Integer> linkedObs = new ArrayList<Integer>(odObs.array.length);

            biasPos = new HashMap<String, Integer>();
            if (odCfg.cfgStations != null)
            {
                String[] stations = odCfg.cfgStations.keySet().toArray(new String[0]);
                for (int i = 0; i < stations.length; i++)
                {
                    for (int j = 0; j < measNames.length; j++)
                    {
                        for (int k = 0; k < odCfg.parameters.size(); k++)
                        {
                            String bias = new StringBuilder(stations[i]).append(measNames[j]).toString();
                            if (bias.equalsIgnoreCase(odCfg.parameters.get(k).name))
                                biasPos.put(bias, k + 6);
                        }
                    }
                }
            }

            int Rsize = 0;
            for (Settings.MeasurementType s: measNames)
                Rsize += odCfg.cfgMeasurements.get(s).error.length;
            Array2DRowRealMatrix R = new Array2DRowRealMatrix(Rsize, Rsize);
            for (int i = 0, j = 0; i < measNames.length; i++)
            {
                Settings.Measurement jm = odCfg.cfgMeasurements.get(measNames[i]);
                for (int k = 0; k < jm.error.length; k++)
                {
                    R.setEntry(j, j, jm.error[k]*jm.error[k]);
                    j++;
                }
            }

            // Create object for each RSO
            ArrayList<ArrayList<SpaceObject>> rsoChoices = new ArrayList<ArrayList<SpaceObject>>(cfgList.size());
            for (Settings cfg: cfgList)
            {
                ArrayList<SpaceObject> rsoList = new ArrayList<SpaceObject>();
                rsoList.add(new SpaceObject(cfg, Rsize));
                rsoChoices.add(rsoList);
            }

            for (int smoothIter = 0; smoothIter < odCfg.estmSmootherIterations; smoothIter++)
            {
                // Re-initialize ICs with previous smoothed results.
                if (smoothIter >= 1)
                {
                    for (Iterator<ArrayList<SpaceObject>> iter = rsoChoices.iterator(); iter.hasNext(); )
                    {
                        SpaceObject currSC = iter.next().get(0);
                        if (currSC.estOutput.size() == 0)
                        {
                            iter.remove();
                            continue;
                        }
                        Estimation.EstimationOutput estOut = currSC.estOutput.get(0);
                        currSC.xhat = new ArrayRealVector(estOut.estimatedState);
                        currSC.xhatPrev = new ArrayRealVector(estOut.estimatedState);
                        if (estOut.estimatedCovariance != null)
                            currSC.P = new Array2DRowRealMatrix(estOut.estimatedCovariance);
                        currSC.dataAssociated = true;
                        tm = estOut.time;
                        currSC.smootherData = new ArrayList<SmootherTimeStep>();
                        currSC.estOutput = new ArrayList<Estimation.EstimationOutput>();
                        currSC.associatedObsIndex = new ArrayList<Integer>();
                        currSC.marginalEvents = new ArrayList<JPDALikelihood>();
                    }
                }

                for (int measIndex = 0, additionalMeas = 0; measIndex < odObs.array.length; measIndex += additionalMeas + 1)
                {
                    additionalMeas = 0;
                    if (linkedObs.contains(measIndex))
                        continue;

                    final AbsoluteDate t0 = tm;
                    tm = odObs.array[measIndex].time;
                    // Determine Number of measurements that need to be collected.
                    while (measIndex + additionalMeas + 1 < odObs.array.length &&
                           odObs.array[measIndex + additionalMeas].equals(odObs.array[measIndex + additionalMeas + 1]))
                        additionalMeas++;

                    // Create new objects with CAR
                    if (activateCAR && rsoChoices.size() == 0)
                    {
                        for (int measNum = 0; measNum <= additionalMeas; measNum++)
                        {
                            ArrayList<Hypothesis> newHypotheses;
                            if (singleObject)
                                newHypotheses = generateOpticalHypotheses(odObs.array[measIndex + measNum], R.getEntry(0, 0),
                                                                          1E-10, R.getEntry(1, 1), 1E-10);
                            else
                                newHypotheses = generateRadarHypotheses(odObs.array[measIndex + measNum], R.getEntry(0, 0), R.getEntry(1, 1),
                                                                        15E-6, 15E-6);

                            ArrayList<SpaceObject> tempHypotheses = new ArrayList<SpaceObject>(newHypotheses.size());
                            for (Hypothesis h: newHypotheses)
                                tempHypotheses.add(new SpaceObject(h.xhat.toArray(), h.P, Rsize, h.weight));
                            rsoChoices.add(tempHypotheses);
                        }
                        continue;
                    }

                    for (ArrayList<SpaceObject> rsoList: rsoChoices)
                    {
                        for (SpaceObject rso: rsoList)
                            rso.marginalEvents.clear();
                    }

                    ArrayList<JPDALikelihood> singleJointEvent = new ArrayList<JPDALikelihood>();
                    ArrayList<ArrayList<JPDALikelihood>> jointEvents = new ArrayList<ArrayList<JPDALikelihood>>();

                    for (int objNum = 0; objNum < rsoChoices.size(); objNum++)
                    {
                        for (Iterator<SpaceObject> iter = rsoChoices.get(objNum).iterator(); iter.hasNext(); )
                        {
                            final SpaceObject currSC = iter.next();
                            final double stepStart = t0.durationFrom(odCfg.propStart);
                            final double stepEnd = tm.durationFrom(odCfg.propStart);

                            if (currSC.dataAssociated)
                            {
                                currSC.dataAssociated = false;
                                currSC.sigma = generateSigmaPoints(currSC.xhat, currSC.P);
                                double[][] sigData = currSC.sigma.getData();
                                for (int j = 6; j < odCfg.parameters.size() + 6; j++)
                                {
                                    Settings.Parameter tempep = odCfg.parameters.get(j - 6);
                                    for (int i = 0; i < numSigmas; i++)
                                        sigData[j][i] = FastMath.min(FastMath.max(sigData[j][i], tempep.min), tempep.max);
                                }
                                currSC.sigma.setSubMatrix(sigData, 0, 0);
                            }
                            else // If data is not associated, keep same sigma points.
                                currSC.sigma.setSubMatrix(currSC.propSigma.getData(), 0, 0);

                            if (FastMath.abs(stepEnd - stepStart) > 1.0E-6)
                            {
                                try
                                {
                                    propagator.propagate(stepStart, currSC.sigma, stepEnd, currSC.propSigma, false);
                                }
                                catch (RuntimeException e)
                                {
                                    // Remove objects that fail to propagate
                                    iter.remove();
                                    continue;
                                }
                            }
                            else
                                currSC.propSigma.setSubMatrix(currSC.sigma.getData(), 0, 0);

                            currSC.xhatPrev = addColumns(currSC.propSigma).mapMultiplyToSelf(weight);
                            currSC.xhat = currSC.xhatPrev;
                            currSC.Pprop = odCfg.getProcessNoiseMatrix(stepEnd - stepStart);
                            for (int i = 0; i < numSigmas; i++)
                            {
                                RealVector y = currSC.propSigma.getColumnVector(i).subtract(currSC.xhatPrev);
                                currSC.Pprop = currSC.Pprop.add(y.outerProduct(y).scalarMultiply(weight));
                            }

                            // Compute predicted measurement, compare to each measurement					    
                            currSC.marginalEvents.add(new JPDALikelihood(objNum, 0, (1 - odCfg.estmDetectionProbability)*currSC.hypothesisWeight));
                            for (int measNum = 0; measNum <= additionalMeas; measNum++)
                            {
                                RealVector rawMeas = updatePrep(currSC, tm, measIndex + measNum);
                                RealVector yhatpre = addColumns(currSC.estimMeas).mapMultiplyToSelf(weight);
                                RealMatrix Pyy = R.copy();
                                for (int i = 0; i < numSigmas; i++)
                                {
                                    RealVector y = currSC.estimMeas.getColumnVector(i).subtract(yhatpre);
                                    Pyy = Pyy.add(y.outerProduct(y).scalarMultiply(weight));
                                }

                                RealVector innov = rawMeas.subtract(yhatpre);
                                RealMatrix mahaTemp = MatrixUtils.createColumnRealMatrix(innov.toArray());
                                RealMatrix maha = mahaTemp.transpose().multiply(MatrixUtils.inverse(Pyy).multiply(mahaTemp));
                                if (Math.sqrt(maha.getEntry(0, 0)) < odCfg.estmGatingThreshold)
                                    currSC.marginalEvents.add(
                                        new JPDALikelihood(objNum, measNum + 1, (1 - new ChiSquaredDistribution(Rsize).cumulativeProbability(
                                                                                     maha.getEntry(0, 0)))*currSC.hypothesisWeight));
                            }
                        }
                    }

                    // JPDA
                    double[][] sumJPDALikelihood = new double[rsoChoices.size()][additionalMeas + 2]; // +1 to account for measurement 0 case
                    for (int objNum = 0; objNum < rsoChoices.size(); objNum++)
                    {
                        for (SpaceObject currSC: rsoChoices.get(objNum))
                        {
                            for (JPDALikelihood like: currSC.marginalEvents)
                                sumJPDALikelihood[objNum][like.measurement] += like.psi;
                        }
                    }

                    jointEvents = JPDAJointEvents(jointEvents, sumJPDALikelihood, singleJointEvent, 0);
                    double[] JPDAProbability = new double[jointEvents.size()];
                    Arrays.fill(JPDAProbability, 1.0);

                    for (int i = 0; i < jointEvents.size(); i++)
                    {
                        for (int j = 0; j < jointEvents.get(i).size(); j++)
                        {
                            double rowSum = 0;
                            for (int measNum = 0; measNum < sumJPDALikelihood[j].length; measNum++)
                                rowSum += sumJPDALikelihood[j][measNum];
                            if (rowSum > 0)
                                JPDAProbability[i] *= jointEvents.get(i).get(j).psi;
                        }	
                    }

                    double JPDAsum = Arrays.stream(JPDAProbability).sum();
                    for (int i = 0; i < jointEvents.size(); i++)
                        JPDAProbability[i] /= JPDAsum;

                    // Identify max probability to see if measurments has been associated
                    int maxProbIndex = 0;
                    for (int i = 1; i < JPDAProbability.length; i++)
                    {
                        if (JPDAProbability[i] > JPDAProbability[maxProbIndex])
                            maxProbIndex = i;
                    }

                    for (int i = 0; i < jointEvents.get(maxProbIndex).size(); i++)
                    {	
                        if (jointEvents.get(maxProbIndex).get(i).measurement != 0)
                        {
                            // Save index
                            int objNum = jointEvents.get(maxProbIndex).get(i).object;
                            int measNum = measIndex + jointEvents.get(maxProbIndex).get(i).measurement - 1;
                            for (int j = 0; j < rsoChoices.get(objNum).size(); j++)
                            {
                                // Indicate all hypotheses for a given object have been associated
                                rsoChoices.get(objNum).get(j).dataAssociated = true;
                                // Associate measurements to hypotheses
                                rsoChoices.get(objNum).get(j).associatedObsIndex.add(measNum);
                            }
                        }
                    }

                    // Combine Hypotheses into Events
                    double[][] betaSatMeas = new double[rsoChoices.size()][additionalMeas + 2]; // +1 to account for measurement 0 case			
                    // Compute beta based on the object & measurement pair
                    for (int i = 0; i < jointEvents.size(); i++)
                    {
                        for (int j = 0; j < jointEvents.get(i).size(); j++)
                            betaSatMeas[jointEvents.get(i).get(j).object][jointEvents.get(i).get(j).measurement] += JPDAProbability[i];
                    }		

                    // Update
                    for (int objNum = 0; objNum < rsoChoices.size(); objNum++)
                    {
                        for (int hypNum = 0; hypNum < rsoChoices.get(objNum).size(); hypNum++)
                        {
                            RealVector betav = new ArrayRealVector(Rsize);
                            RealMatrix betavvt = new Array2DRowRealMatrix(Rsize, Rsize);
                            SpaceObject currSC = rsoChoices.get(objNum).get(hypNum);

                            // Update hypothesis weight
                            if (rsoChoices.get(objNum).size() > 1)
                            {
                                double hypWeight = 0;			
                                for (int counter = 0; counter < currSC.marginalEvents.size(); counter++)
                                {
                                    int measNum = currSC.marginalEvents.get(counter).measurement;
                                    hypWeight += betaSatMeas[objNum][measNum]/sumJPDALikelihood[objNum][measNum]*
                                        currSC.marginalEvents.get(counter).psi;
                                }
                                currSC.hypothesisWeight = hypWeight;
                            }

                            for (int i = 0; i < currSC.marginalEvents.size(); i++)
                            {
                                if (currSC.marginalEvents.get(i).measurement > 0)
                                {
                                    double betaTemp = 0;
                                    RealVector rawMeas = updatePrep(currSC, tm, measIndex + currSC.marginalEvents.get(i).measurement - 1);
                                    RealVector yhatpre = addColumns(currSC.estimMeas).mapMultiplyToSelf(weight);
                                    RealVector innov = quadCheck(rawMeas, yhatpre);
                                    if (currSC.marginalEvents.get(i).psi != 0)
                                    {
                                        betaTemp = betaSatMeas[objNum][currSC.marginalEvents.get(i).measurement]/
                                            sumJPDALikelihood[objNum][currSC.marginalEvents.get(i).measurement]*
                                            currSC.marginalEvents.get(i).psi/currSC.hypothesisWeight;
                                    }
                                    betav = betav.add(innov.mapMultiply(betaTemp));
                                    betavvt = betavvt.add(innov.outerProduct(innov).scalarMultiply(betaTemp));
                                }
                            }

                            RealMatrix Pyy = R.copy();
                            RealMatrix Pxy = new Array2DRowRealMatrix(numStates, Rsize);
                            RealVector yhatpre = addColumns(currSC.estimMeas).mapMultiplyToSelf(weight);
                            for (int i = 0; i < numSigmas; i++)
                            {
                                RealVector y = currSC.estimMeas.getColumnVector(i).subtract(yhatpre);
                                Pyy = Pyy.add(y.outerProduct(y).scalarMultiply(weight));
                                Pxy = Pxy.add(currSC.propSigma.getColumnVector(i).subtract(currSC.xhatPrev).outerProduct(y).scalarMultiply(weight));
                            }

                            RealMatrix K = Pxy.multiply(MatrixUtils.inverse(Pyy));
                            RealMatrix Ptilda = K.multiply((betavvt.subtract(betav.outerProduct(betav))).multiply(K.transpose()));
                            currSC.xhat = currSC.xhatPrev.add(odCfg.parameterMatrix.multiply(K).operate(betav));
                            currSC.P = currSC.Pprop.subtract(odCfg.parameterMatrix.multiply(K.multiply(Pyy.multiply(K.transpose()))).subtract(Ptilda));

                            // For estimated parameters, return xhat values to max/min bounds. If far outside the bounds, can result in
                            // singular matrix when taking the inverse for smoothing.
                            for (int j = 6; j < odCfg.parameters.size() + 6; j++)
                            {
                                Settings.Parameter tempep = odCfg.parameters.get(j - 6);
                                currSC.xhat.setEntry(j, FastMath.min(FastMath.max(currSC.xhat.getEntry(j), tempep.min), tempep.max));
                            }

                            // If data associated, save data and compute smoother values.
                            if (currSC.dataAssociated)
                            {
                                double[] pv = currSC.xhat.toArray();
                                CartesianOrbit orb = new CartesianOrbit(new PVCoordinates(new Vector3D(pv[0],pv[1],pv[2]),
                                                                                          new Vector3D(pv[3],pv[4],pv[5])),
                                                                        odCfg.propInertialFrame, tm, Constants.EGM96_EARTH_MU);
                                ssta[0] = new SpacecraftState(orb, odCfg.rsoMass);

                                Estimation.EstimationOutput odout = new Estimation.EstimationOutput(odObs.array[measIndex].time,
                                                                                                    odObs.array[measIndex].station);
                                currSC.estOutput.add(odout);
                                odout.estimatedState = pv;
                                odout.propagatedCovariance = getLowerTriangle(currSC.Pprop);
                                odout.innovationCovariance = getLowerTriangle(Pyy);
                                odout.estimatedCovariance = getLowerTriangle(currSC.P);
                                odout.preFit = yhatpre.toArray();
                                if (singleObject)
                                    odout.postFit = odObs.array[measIndex].helpers[0].estimate(0, 0, ssta).getEstimatedValue();
                                else
                                    odout.postFit = new double[]{odObs.array[measIndex].helpers[0].estimate(0, 0, ssta).getEstimatedValue()[0],
                                        odObs.array[measIndex].helpers[1].estimate(0, 0, ssta).getEstimatedValue()[0]};

                                // Generate post sigma points
                                Array2DRowRealMatrix postSigma = generateSigmaPoints(currSC.xhat, currSC.P);
                                double[][] sigData = postSigma.getData();
                                for (int j = 6; j < odCfg.parameters.size() + 6; j++)
                                {
                                    Settings.Parameter tempep = odCfg.parameters.get(j - 6);
                                    for (int i = 0; i < numSigmas; i++)
                                        sigData[j][i] = FastMath.min(FastMath.max(sigData[j][i], tempep.min), tempep.max);
                                }
                                postSigma.setSubMatrix(sigData, 0, 0);

                                // Store Smoother Data
                                SmootherTimeStep smout = new SmootherTimeStep();
                                smout.xpre = MatrixUtils.createColumnRealMatrix(currSC.xhatPrev.toArray());
                                smout.xpost = MatrixUtils.createColumnRealMatrix(currSC.xhat.toArray());
                                smout.Ppre = currSC.Pprop;
                                smout.Ppost = currSC.P;
                                smout.sigPre = MatrixUtils.createRealMatrix(currSC.propSigma.getData());
                                smout.sigPost = MatrixUtils.createRealMatrix(postSigma.getData());
                                smout.xstar = smout.xpost;
                                smout.Pstar = smout.Ppost;
                                smout.tmSmoother = tm;
                                smout.measObjs = odObs.array[measIndex].helpers;
                                currSC.smootherData.add(smout);
                            }
                        }
                    }

                    // Merge hypotheses and normalize weights
                    for (ArrayList<SpaceObject> rsoList: rsoChoices)
                    {
                        double sumWeights = 0;
                        mergeHypotheses(rsoList);
                        for (SpaceObject rso: rsoList)
                            sumWeights += rso.hypothesisWeight;
                        for (SpaceObject rso: rsoList)
                            rso.hypothesisWeight /= sumWeights;
                    }
                }

                // Smooth data		
                for (ArrayList<SpaceObject> rsoList: rsoChoices)
                {
                    for (int hypNum = 0; hypNum < rsoList.size(); hypNum++)
                    {
                        SpaceObject currSC = rsoList.get(hypNum);
                        currSC.isConsistent = true;
                        int smSize = currSC.smootherData.size() - 1;

                        for (int i = 0; i < smSize; i++)
                        {
                            SmootherTimeStep smDatak1 = currSC.smootherData.get(smSize - i), smDatak = currSC.smootherData.get(smSize - i - 1);
                            RealMatrix Csmoother = new Array2DRowRealMatrix(numStates, numStates);
                            RealMatrix Asmoother = new Array2DRowRealMatrix(numStates, numStates);
                            for (int j = 0; j < numSigmas; j++) 
                                Csmoother = Csmoother.add(smDatak.sigPost.getColumnMatrix(j).subtract(smDatak.xpost).multiply(
                                                              smDatak1.sigPre.getColumnMatrix(j).subtract(smDatak1.xpre)
                                                              .transpose()).scalarMultiply(weight));
                            Asmoother = Csmoother.multiply(MatrixUtils.inverse(smDatak1.Ppre));
                            smDatak.xstar = smDatak.xpost.add(Asmoother.multiply(smDatak1.xstar.subtract(smDatak1.xpre)));
                            smDatak.Pstar = smDatak.Ppost.add(Asmoother.multiply((smDatak1.Pstar.subtract(smDatak1.Ppre))
                                                                                 .multiply(Asmoother.transpose())));
                            currSC.smootherData.set(smSize - i - 1,smDatak);
                        }

                        // Compute McReynolds consistency or merge all hypotheses into one hypothesis
                        if (!rsoList.get(0).fromCAR)
                        {
                            for (int i = 0; i < currSC.smootherData.size()-1;i++)
                            {
                                SmootherTimeStep smDatak = currSC.smootherData.get(smSize - i - 1);
                                RealMatrix delx = smDatak.xpost.subtract(smDatak.xstar), delP = smDatak.Ppost.subtract(smDatak.Pstar);
                                for (int j = 0; j < 5; j++)
                                {			
                                    if (Math.abs(delx.getEntry(j, 0))/Math.sqrt(Math.abs(delP.getEntry(j, j))) >= 3)
                                        currSC.isConsistent = false;
                                }
                            }

                            SpacecraftState[] sms = new SpacecraftState[1];
                            for (int j = 0; j < currSC.smootherData.size();j++)
                            {
                                double[] pv = currSC.smootherData.get(j).xstar.getColumn(0);
                                tm = currSC.smootherData.get(j).tmSmoother;
                                sms[0] = new SpacecraftState(new CartesianOrbit(new PVCoordinates(new Vector3D(pv[0], pv[1], pv[2]),
                                                                                                  new Vector3D(pv[3], pv[4], pv[5])),
                                                                                odCfg.propInertialFrame, tm, Constants.EGM96_EARTH_MU));
                                currSC.estOutput.get(j).estimatedState = currSC.smootherData.get(j).xstar.getColumn(0);
                                currSC.estOutput.get(j).estimatedCovariance = getLowerTriangle(currSC.smootherData.get(j).Pstar);
                                // Compute smoothed residuals
                                if (singleObject)
                                    currSC.estOutput.get(j).postFit = currSC.smootherData.get(j).measObjs[0].estimate(0, 0, sms).getEstimatedValue();
                                else
                                {
                                    currSC.estOutput.get(j).postFit = new double[]{currSC.smootherData.get(j).measObjs[0].estimate(0, 0, sms).
                                        getEstimatedValue()[0],	currSC.smootherData.get(j).measObjs[1].estimate(0, 0, sms).getEstimatedValue()[0]};
                                }
                            }
                        }
                        else
                        {
                            rsoList.get(0).fromCAR = false;
                            currSC.isConsistent = false;
                            break;
                        }
                    }
                }

                // Check McReynolds and potentially break out
                boolean trackPromoted = false;
                for (Iterator<ArrayList<SpaceObject>> iter = rsoChoices.iterator(); iter.hasNext(); )
                {
                    SpaceObject currSC = iter.next().get(0);
                    if (currSC.isConsistent)
                    {
                        trackPromoted = true;
                        for (int idx = 0; idx < currSC.associatedObsIndex.size(); idx++)
                        {
                            RealMatrix Ptemp = new Array2DRowRealMatrix(currSC.estOutput.get(idx).estimatedCovariance);
                            currSC.xhatPrev = new ArrayRealVector(currSC.estOutput.get(idx).estimatedState);
                            currSC.propSigma = generateSigmaPoints(currSC.xhatPrev, Ptemp);

                            RealMatrix Pyy = R.copy();
                            RealVector yhatpre = addColumns(currSC.estimMeas).mapMultiplyToSelf(weight);
                            RealVector rawMeas = updatePrep(currSC, tm, currSC.associatedObsIndex.get(idx));
                            for (int i = 0; i < numSigmas; i++)
                            {
                                RealVector y = currSC.estimMeas.getColumnVector(i).subtract(yhatpre);
                                Pyy = Pyy.add(y.outerProduct(y).scalarMultiply(weight));
                            }
                            currSC.estOutput.get(idx).innovationCovariance = getLowerTriangle(Pyy);

                            // Add each measurement to the list of associated measurements.
                            linkedObs.add(currSC.associatedObsIndex.get(idx));
                        }

                        promotedTracks.add(currSC);
                        iter.remove();
                    }
                }

                if (trackPromoted)
                {
                    if (odObs.array.length - linkedObs.size() <= 5)
                        break;
                    if (rsoChoices.size() == 0)
                    {
                        smoothIter = -1;
                        activateCAR = true;
                    }
                }
            }

            // Save output
            multiOutput.estOutput = new ArrayList<ArrayList<Estimation.EstimationOutput>>(promotedTracks.size());
            multiOutput.associatedObs = new ArrayList<ArrayList<Integer>>(odObs.array.length);
            multiOutput.unassociatedObs = new ArrayList<Integer>(odObs.array.length);
            for (SpaceObject obj: promotedTracks)
            {
                multiOutput.estOutput.add(obj.estOutput);
                multiOutput.associatedObs.add(obj.associatedObsIndex);
            }

            for (int idx = 0; idx < odObs.array.length; idx++)
            {
                if (!linkedObs.contains(idx))
                    multiOutput.unassociatedObs.add(idx);
            }
        }

        private	ArrayList<Hypothesis> generateOpticalHypotheses(Measurements.Measurement obs, double sigmaRA, double sigmaRAd,
                                                                double sigmaDec, double sigmaDecd)
        {
            double[] x0 = odCfg.getInitialState();
            double RA = obs.values[0], RA_d = obs.angleRates[0], Dec = obs.values[1], Dec_d = obs.angleRates[1];
            TimeStampedPVCoordinates station = odCfg.stations.get(obs.station).getBaseFrame().getPVCoordinates(obs.time, odCfg.propInertialFrame);
            RealVector meanTemp = new ArrayRealVector(numStates);
            RealMatrix covarTemp = new DiagonalMatrix(numStates);
            ArrayList<Hypothesis> objectHypotheses = new ArrayList<Hypothesis>();

            ArrayList <CAR.CARGaussianElement> CARGaussians = new CAR(RA, Dec, RA_d, Dec_d, station, 100000.0, 100.0, 10000.0,
                                                                      17000000, 37000000, 0.05, singleObject).getCAR();
            for (int i = 0; i < CARGaussians.size(); i++)
            {
                if (CARGaussians.get(i).ordinateStd < 0)
                    continue;
                meanTemp.setEntry(0, CARGaussians.get(i).abscissaMean);	meanTemp.setEntry(1, RA); meanTemp.setEntry(2, Dec);
                meanTemp.setEntry(3, CARGaussians.get(i).ordinateMean);	meanTemp.setEntry(4, RA_d); meanTemp.setEntry(5, Dec_d);
                covarTemp.setEntry(0, 0, Math.pow(CARGaussians.get(i).abscissaStd, 2));	covarTemp.setEntry(1, 1, sigmaRA);
                covarTemp.setEntry(2, 2, sigmaDec); covarTemp.setEntry(3, 3, Math.pow(CARGaussians.get(i).ordinateStd, 2));
                covarTemp.setEntry(4, 4, sigmaRAd); covarTemp.setEntry(5, 5, sigmaDecd);
                if (numStates >= 7)
                {
                    meanTemp.setEntry(6, x0[6]);
                    covarTemp.setEntry(6, 6, Math.pow(odCfg.estmCovariance[6], 2));
                    if (numStates >= 8)
                    {
                        meanTemp.setEntry(7, x0[7]);
                        covarTemp.setEntry(7, 7, Math.pow(odCfg.estmCovariance[7], 2));
                    }
                }

                Array2DRowRealMatrix sigma = generateSigmaPoints(meanTemp, covarTemp);
                for (int j = 0; j < numSigmas; j++)
                    sigma.setColumn(j, rangeRaDecToCart(sigma.getColumnVector(j), obs.time, obs.station));

                RealMatrix Pxx = new Array2DRowRealMatrix(numStates, numStates);
                RealVector xhat = addColumns(sigma).mapMultiplyToSelf(0.5/numStates);
                for (int j = 0; j < numSigmas; j++)
                {
                    RealVector xhatdiff = sigma.getColumnVector(j).subtract(xhat);
                    Pxx = Pxx.add(xhatdiff.outerProduct(xhatdiff).scalarMultiply(0.5/numStates));
                }
                objectHypotheses.add(new Hypothesis(xhat, Pxx, CARGaussians.get(i).weight));
            }
            return(objectHypotheses);
        }

        private ArrayList<Hypothesis> generateRadarHypotheses(Measurements.Measurement obs, double sigmaRange, double sigmaRR,
                                                              double sigmaRA, double sigmaDec)
        {
            double[] x0 = odCfg.getInitialState();
            double range = obs.values[0], rangeRate = obs.values[1], ra = obs.angleRates[0], dec = obs.angleRates[1];
            TimeStampedPVCoordinates station = odCfg.stations.get(obs.station).getBaseFrame().getPVCoordinates(obs.time, odCfg.propInertialFrame);
            RealVector meanTemp = new ArrayRealVector(numStates);
            RealMatrix covarTemp = new DiagonalMatrix(numStates);
            ArrayList<Hypothesis> objectHypotheses = new ArrayList<Hypothesis>();

            ArrayList <CAR.CARGaussianElement> CARGaussians = new CAR(ra, dec, range, rangeRate, station, 5E-6, 5E-6, 5E-6,
                                                                      17000000, 37000000, 0.1, singleObject).getCAR();
            for (int i = 0; i < CARGaussians.size(); i++)
            {
                if (CARGaussians.get(i).ordinateStd < 0)
                    continue;
                meanTemp.setEntry(0, range); meanTemp.setEntry(1, ra); meanTemp.setEntry(2, dec); meanTemp.setEntry(3, rangeRate);
                meanTemp.setEntry(4, CARGaussians.get(i).abscissaMean); meanTemp.setEntry(5, CARGaussians.get(i).ordinateMean);
                covarTemp.setEntry(0, 0, sigmaRange); covarTemp.setEntry(1, 1, sigmaRA); covarTemp.setEntry(2, 2, sigmaDec);
                covarTemp.setEntry(3, 3, sigmaRR); covarTemp.setEntry(4, 4, Math.pow(CARGaussians.get(i).abscissaStd, 2));
                covarTemp.setEntry(5, 5, Math.pow(CARGaussians.get(i).ordinateStd, 2));
                if (numStates >= 7)
                {
                    meanTemp.setEntry(6, x0[6]);
                    covarTemp.setEntry(6, 6, Math.pow(odCfg.estmCovariance[6], 2));
                    if (numStates >= 8)
                    {
                        meanTemp.setEntry(7, x0[7]);
                        covarTemp.setEntry(7, 7, Math.pow(odCfg.estmCovariance[7], 2));
                    }
                }

                Array2DRowRealMatrix sigma = generateSigmaPoints(meanTemp, covarTemp);
                for (int j = 0; j < numSigmas; j++)
                    sigma.setColumn(j, rangeRaDecToCart(sigma.getColumnVector(j), obs.time, obs.station));

                RealMatrix Pxx = new Array2DRowRealMatrix(numStates, numStates);
                RealVector xhat = addColumns(sigma).mapMultiplyToSelf(0.5/numStates);
                for (int j = 0; j < numSigmas; j++)
                {
                    RealVector xhatdiff = sigma.getColumnVector(j).subtract(xhat);
                    Pxx = Pxx.add(xhatdiff.outerProduct(xhatdiff).scalarMultiply(0.5/numStates));
                }
                objectHypotheses.add(new Hypothesis(xhat, Pxx, CARGaussians.get(i).weight));
            }
            return(objectHypotheses);
        }

        private void mergeHypotheses(ArrayList<SpaceObject> hypotheses)
        {
            double maxWeight = 0.0;
            for (SpaceObject rso: hypotheses)
            {
                if (rso.hypothesisWeight > maxWeight)
                    maxWeight = rso.hypothesisWeight;
            }
            maxWeight *= 1E-3;

            for (int hyp1 = 0; hyp1 < hypotheses.size(); hyp1++)
            {
                SpaceObject obj1 = hypotheses.get(hyp1);
                for (Iterator<SpaceObject> iter2 = hypotheses.iterator(); iter2.hasNext(); )
                {
                    SpaceObject obj2 = iter2.next();
                    if (obj1 == obj2)
                        continue;
                    if (obj2.hypothesisWeight < maxWeight)
                    {
                        iter2.remove();
                        continue;
                    }
                    RealMatrix mahaTemp = MatrixUtils.createColumnRealMatrix(obj1.xhat.subtract(obj2.xhat).toArray());
                    RealMatrix maha1 = mahaTemp.transpose().multiply(MatrixUtils.inverse(obj1.P).multiply(mahaTemp));
                    RealMatrix maha2 = mahaTemp.transpose().multiply(MatrixUtils.inverse(obj2.P).multiply(mahaTemp));

                    // If Mahalanobis distance low, combine distributions.
                    if (Math.min(maha1.getEntry(0, 0), maha2.getEntry(0, 0)) < 1)
                    {
                        RealVector xhatTemp =  obj1.xhat.mapMultiply(obj1.hypothesisWeight).add(obj2.xhat.mapMultiply(obj2.hypothesisWeight))
                            .mapMultiply(1/(obj1.hypothesisWeight + obj2.hypothesisWeight));
                        RealMatrix matrixTemp = MatrixUtils.createColumnRealMatrix(obj1.xhat.toArray());
                        RealMatrix PTemp = obj1.P.add(matrixTemp.multiply(matrixTemp.transpose()))
                            .scalarMultiply(obj1.hypothesisWeight/(obj1.hypothesisWeight + obj2.hypothesisWeight));
                        matrixTemp = MatrixUtils.createColumnRealMatrix(obj2.xhat.toArray());
                        PTemp = PTemp.add(obj2.P.add(matrixTemp.multiply(matrixTemp.transpose()))
                                          .scalarMultiply(obj2.hypothesisWeight/(obj1.hypothesisWeight + obj2.hypothesisWeight)));
                        matrixTemp = MatrixUtils.createColumnRealMatrix(xhatTemp.toArray());
                        obj1.xhat = xhatTemp;
                        obj1.P = PTemp.subtract(matrixTemp.multiply(matrixTemp.transpose()));
                        obj1.hypothesisWeight += obj2.hypothesisWeight;
                        iter2.remove();
                    }
                }
            }
        }

        private Array2DRowRealMatrix generateSigmaPoints(RealVector mean, RealMatrix cov)
        {
            RealMatrix Ptemp;
            if (cov.getRowDimension() == cov.getColumnDimension())
                Ptemp = cov.scalarMultiply(numStates);
            else
            {
                Ptemp = new Array2DRowRealMatrix(numStates, numStates);
                for (int i = 0, k = 0; i < numStates; i++)
                {
                    for (int j = 0; j < i + 1; j++, k++)
                    {
                        double value = cov.getEntry(k, 0)*numStates;
                        Ptemp.setEntry(i, j, value);
                        Ptemp.setEntry(j, i, value);
                    }
                }
            }

            RealMatrix sqrP = new CholeskyDecomposition(Ptemp.add(Ptemp.transpose()).scalarMultiply(0.5), 1E-6, 1E-16).getL();
            Array2DRowRealMatrix sigma = new Array2DRowRealMatrix(numStates, numSigmas);
            for (int j = 0; j < numStates; j++)
            {
                sigma.setColumnVector(j, mean.add(sqrP.getColumnVector(j)));
                sigma.setColumnVector(numStates + j, mean.subtract(sqrP.getColumnVector(j)));
            }
            return(sigma);
        }

        private	double[] rangeRaDecToCart(RealVector rRaDec, AbsoluteDate date, String stat)
        {
            double Range = rRaDec.getEntry(0), Range_d = rRaDec.getEntry(3);
            double RA = rRaDec.getEntry(1), RA_d = rRaDec.getEntry(4);
            double Dec = rRaDec.getEntry(2), Dec_d = rRaDec.getEntry(5);

            // Compute inertial pos/vel relative to station
            Vector3D topoPos = new Vector3D(new double[] {Range*Math.cos(Dec)*Math.cos(RA), Range*Math.cos(Dec)*Math.sin(RA), Range*Math.sin(Dec)});
            Vector3D topoVel = new Vector3D(new double[] {Range_d*Math.cos(Dec)*Math.cos(RA) - Range*Math.sin(Dec)*Math.cos(RA)*Dec_d -
                                                          Range*Math.cos(Dec)*Math.sin(RA)*RA_d, Range_d*Math.cos(Dec)*Math.sin(RA) -
                                                          Range*Math.sin(Dec)*Math.sin(RA)*Dec_d +
                                                          Range*Math.cos(Dec)*Math.cos(RA)*RA_d, Range_d*Math.sin(Dec) + Range*Math.cos(Dec)*Dec_d});

            // Add station coords to observation coords.
            Vector3D pos = new Vector3D(1, topoPos, 1, odCfg.stations.get(stat).getBaseFrame().
                                        getPVCoordinates(date, odCfg.propInertialFrame).getPosition());
            Vector3D vel = new Vector3D(1, topoVel, 1, odCfg.stations.get(stat).getBaseFrame().
                                        getPVCoordinates(date, odCfg.propInertialFrame).getVelocity());

            if (rRaDec.getDimension() >= 8)
                return(new double[] {pos.getX(), pos.getY(), pos.getZ(), vel.getX(), vel.getY(), vel.getZ(),
                                     rRaDec.getEntry(6), rRaDec.getEntry(7)});
            if (rRaDec.getDimension() >= 7)
                return(new double[] {pos.getX(), pos.getY(), pos.getZ(), vel.getX(), vel.getY(), vel.getZ(), rRaDec.getEntry(6)});
            return(new double[] {pos.getX(), pos.getY(), pos.getZ(), vel.getX(), vel.getY(), vel.getZ()});
        }

        private RealVector quadCheck(RealVector measurement, RealVector state)
        {
            // Innovations QuadCheck
            RealVector innov = measurement.subtract(state);
            if (singleObject && innov.getEntry(0) > Math.PI)
                innov.setEntry(0, innov.getEntry(0) - 2*Math.PI);
            else if (singleObject && innov.getEntry(0) < -Math.PI)
                innov.setEntry(0, innov.getEntry(0) + 2*Math.PI);
            return(innov);
        }

        private ArrayList<ArrayList<JPDALikelihood>> JPDAJointEvents(ArrayList<ArrayList<JPDALikelihood>> jointEvents, double[][] marginalEvents,
                                                                     ArrayList<JPDALikelihood> singleJointEvent, int objNum)
        {
            for (int measNum = 0; measNum < marginalEvents[objNum].length; measNum++)
            {	
                // Copy Single Joint Event into Single Joint Event Temp
                ArrayList<JPDALikelihood> singleJointEventTemp = new ArrayList<JPDALikelihood>();
                for (int m = 0; m < singleJointEvent.size(); m++)
                    singleJointEventTemp.add(singleJointEvent.get(m));

                // Add event if that measurement has not been used
                boolean skip = false;
                for (int m = 0; m < singleJointEvent.size(); m++)
                {
                    if (measNum == singleJointEventTemp.get(m).measurement && measNum != 0)
                        skip = true;
                }
                if (skip)
                    continue;

                singleJointEventTemp.add(new JPDALikelihood(objNum, measNum, marginalEvents[objNum][measNum]));
                if (marginalEvents.length == objNum + 1)
                    jointEvents.add(singleJointEventTemp);
                else
                    jointEvents = JPDAJointEvents(jointEvents, marginalEvents, singleJointEventTemp, objNum + 1);
            }
            return(jointEvents);
        }

        private RealVector updatePrep(SpaceObject currSC, AbsoluteDate tm, int measIndex)
        {		
            // Transforms measurement Sig Points and stores in currSC. Also returns rawMeas
            RealVector rawMeas = null;
            for (int i = 0; i < numSigmas; i++)
            {
                double[] pv = currSC.propSigma.getColumn(i);
                final SpacecraftState[] ssta = new SpacecraftState[1];
                CartesianOrbit orb = new CartesianOrbit(new PVCoordinates(new Vector3D(pv[0], pv[1], pv[2]), new Vector3D(pv[3], pv[4], pv[5])),
                                                        odCfg.propInertialFrame, tm, Constants.EGM96_EARTH_MU);
                ssta[0] = new SpacecraftState(orb, odCfg.rsoMass);

                if (singleObject)
                {
                    currSC.estimMeas.setColumn(i, odObs.array[measIndex].helpers[0].estimate(0, 0, ssta).getEstimatedValue());
                    if (rawMeas == null)
                        rawMeas = new ArrayRealVector(odObs.array[measIndex].helpers[0].getObservedValue());
                }
                else
                {
                    currSC.estimMeas.setEntry(0, i, odObs.array[measIndex].helpers[0].estimate(0, 0, ssta).getEstimatedValue()[0]);
                    currSC.estimMeas.setEntry(1, i, odObs.array[measIndex].helpers[1].estimate(0, 0, ssta).getEstimatedValue()[0]);
                    if (rawMeas == null)
                        rawMeas = new ArrayRealVector(new double[]{odObs.array[measIndex].helpers[0].getObservedValue()[0],
                                                                   odObs.array[measIndex].helpers[1].getObservedValue()[0]});
                }
            }

            // Does not enter if meas is pos/vel. Otherwise will enter.
            if (odObs.array[measIndex].station != null)
            {
                String name = new StringBuilder(odObs.array[measIndex].station).append(measNames[0]).toString();
                Integer pos = biasPos.get(name);
                if (pos != null)
                {
                    if (measNames.length == 2)
                    {
                        name = new StringBuilder(odObs.array[measIndex].station).append(measNames[1]).toString();
                        rawMeas = rawMeas.subtract(new ArrayRealVector(new double[]{currSC.xhatPrev.getEntry(pos),
                                                                                    currSC.xhatPrev.getEntry(pos)}));
                    }
                    else
                        rawMeas = rawMeas.subtract(new ArrayRealVector(new double[]{currSC.xhatPrev.getEntry(pos)}));
                }
            }

            return(rawMeas);
        }

        private ArrayRealVector addColumns(RealMatrix mat)
        {
            double[][] arr = mat.getData();
            int m = mat.getRowDimension();
            int n = mat.getColumnDimension();
            ArrayRealVector out = new ArrayRealVector(m);
            for (int j = 0; j < m; j++)
            {
                double sum = 0.0;
                for (int i = 0; i < n; i++)
                    sum += arr[j][i];
                out.setEntry(j, sum);
            }
            return(out);
        }

        private double[] getLowerTriangle(RealMatrix mat)
        {
            int m = mat.getRowDimension();
            double[] out = new double[(int)(m*(m + 1)/2)];
            for (int i = 0, k = 0; i < m; i++)
            {
                for (int j = 0; j <= i; j++)
                    out[k++] = mat.getEntry(i, j);
            }
            return(out);
        }
    }
}
