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
	    ObservedMeasurement measObjsSmoother;
	}

	private class JPDALikelihoods
	{
	    double psi;
	    int object;
	    int measurement;
	}

	private class OneObject
	{
	    double[] xInitial;
	    RealMatrix PInitial;
	    RealMatrix P;
	    RealMatrix psdCorr;
	    Array2DRowRealMatrix sigma;
	    Array2DRowRealMatrix propSigma;
	    Array2DRowRealMatrix estimMeas;

	    RealMatrix Pprop;
	    ArrayRealVector xhat;
	    RealVector xhatPrev;
	    double hypothesisWeight;

	    ArrayList<SmootherTimeStep> smootherData;
	    ArrayList<JPDALikelihoods> marginalEvents;
	    ArrayList<Integer> associatedObsIndex;
	    ArrayList<Estimation.EstimationOutput> estOutput;

	    boolean dataAssociated;
	    boolean McReynoldsConsistencyPass;
	    boolean fromCAR;

	    public OneObject(Settings Config, int numStates, int numSigmas, int Rsize)
	    {
		xInitial = Config.getInitialState();
		PInitial = new DiagonalMatrix(Config.estmCovariance);
		P = PInitial;
		psdCorr = MatrixUtils.createRealIdentityMatrix(P.getRowDimension()).scalarMultiply(1.0E-6);
		sigma = new Array2DRowRealMatrix(numStates, numSigmas);
		propSigma = new Array2DRowRealMatrix(numStates, numSigmas);
		estimMeas = new Array2DRowRealMatrix(Rsize, numSigmas);
		Pprop = null;
		xhat = new ArrayRealVector(xInitial);
		xhatPrev = new ArrayRealVector(xInitial);
		smootherData = new ArrayList<SmootherTimeStep>();
		marginalEvents = new ArrayList<JPDALikelihoods>();
		associatedObsIndex = new ArrayList<Integer>();
		estOutput = new ArrayList<Estimation.EstimationOutput>();
		dataAssociated = true;
		hypothesisWeight = 1;
		fromCAR = false;
	    }

	    public OneObject(double[] x, RealMatrix Covar, int numStates, int numSigmas, int Rsize, double hypWeight)
	    {
		xInitial = x;
		PInitial = Covar;
		P = PInitial;
		psdCorr = MatrixUtils.createRealIdentityMatrix(P.getRowDimension()).scalarMultiply(1.0E-6);
		sigma = new Array2DRowRealMatrix(numStates, numSigmas);
		propSigma = new Array2DRowRealMatrix(numStates, numSigmas);
		estimMeas = new Array2DRowRealMatrix(Rsize, numSigmas);
		Pprop = null;
		xhat = new ArrayRealVector(xInitial);
		xhatPrev = new ArrayRealVector(xInitial);
		smootherData = new ArrayList<SmootherTimeStep>();
		marginalEvents = new ArrayList<JPDALikelihoods>();
		associatedObsIndex = new ArrayList<Integer>();
		estOutput = new ArrayList<Estimation.EstimationOutput>();
		dataAssociated = true;
		hypothesisWeight = hypWeight;
		fromCAR = true;
	    }
	}

	private class Hypothesis
	{
	    ArrayRealVector xhat;
 	    RealMatrix P;    
 	    double weight;
	}

	private void determineOrbit()
	{
	    boolean activateCAR = false;
	    final int numStates = odCfg.parameters.size() + 6;
	    final int numSigmas = 2*numStates;
	    final double weight = 0.5/numStates;
	    AbsoluteDate tm = odCfg.propStart;
	    final SpacecraftState[] ssta = new SpacecraftState[1];
	    final ManualPropagation propagator = new ManualPropagation(odCfg);
	    ArrayList<OneObject> promotedTracks = new ArrayList<OneObject>();
	    final ArrayList<Integer> removeObs = new ArrayList<Integer>();

	    HashMap<String, Integer> biasPos = new HashMap<String, Integer>();
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
	    ArrayList<ArrayList<OneObject>> rsoChoices = new ArrayList<ArrayList<OneObject>>(cfgList.size());
	    for (int objNum = 0; objNum < cfgList.size(); objNum++)
	    {
		rsoChoices.add(new ArrayList<OneObject>());
		rsoChoices.get(objNum).add(new OneObject(cfgList.get(objNum), numStates, numSigmas, Rsize));
	    }

	    for (int smoothIter = 0; smoothIter < odCfg.estmSmootherIterations; smoothIter++)
	    {
		// Re-initialize ICs with previous smoothed results.
		if (smoothIter >= 1)
		{
		    for (int objNum = 0; objNum < rsoChoices.size(); objNum++)
		    {
			// TODO Hypotheses need to be merged by this point, since it pulls the 0th hypothesis
			OneObject currSC = rsoChoices.get(objNum).get(0);
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
			currSC.marginalEvents = new ArrayList<JPDALikelihoods>();
		    }
		}

		int measIndex = -1, additionalMeas = 0;
		while (true)
		{
		    final AbsoluteDate t0 = tm;
		    measIndex += additionalMeas + 1;
		    additionalMeas = 0;
		    if (removeObs.contains(measIndex))
			continue;

		    if (measIndex < odObs.array.length)
		    {
			tm = odObs.array[measIndex].time;
			// Determine Number of measurements that need to be collected.
			while (true)
			{
			    if (measIndex + additionalMeas + 1 < odObs.array.length &&
				odObs.array[measIndex + additionalMeas].equals(odObs.array[measIndex + additionalMeas + 1]))
				additionalMeas++;
			    else
				break;
			}
		    }
		    else
			break;

		    // Create new objects with CAR
		    if (activateCAR && rsoChoices.size() == 0)
		    {
			for (int measNum = 0; measNum <= additionalMeas; measNum++)
			{
			    ArrayList<Hypothesis> newHypotheses = new ArrayList<Hypothesis>();
			    if (singleObject)
				newHypotheses = generateOpticalHypotheses(odObs.array[measIndex + measNum], R.getEntry(0, 0), 1E-10,
									  R.getEntry(1, 1), 1E-10, numStates, numSigmas);
			    else
				newHypotheses = generateRadarHypotheses(odObs.array[measIndex + measNum], R.getEntry(0, 0), R.getEntry(1, 1),
									15E-6, 15E-6, numStates, numSigmas);

			    ArrayList<OneObject> tempHypotheses = new ArrayList<OneObject>(newHypotheses.size());
			    for (int hypNum = 0; hypNum < newHypotheses.size(); hypNum++)
			    {
				Hypothesis temp = newHypotheses.get(hypNum);
				OneObject tempState = new OneObject(temp.xhat.toArray(), temp.P, numStates, numSigmas, Rsize, temp.weight);
				tempHypotheses.add(tempState);
			    }
			    rsoChoices.add(tempHypotheses);
			}
			continue;
		    }

		    for (int objNum = 0; objNum < rsoChoices.size(); objNum++)
		    {
			for (int hypNum = 0; hypNum < rsoChoices.get(objNum).size(); hypNum++)
			    rsoChoices.get(objNum).get(hypNum).marginalEvents.clear();
		    }

		    ArrayList<ArrayList<JPDALikelihoods>> jointEvents = new ArrayList<ArrayList<JPDALikelihoods>>();
		    ArrayList<JPDALikelihoods> singleJointEvent = new ArrayList<JPDALikelihoods>();
		    ArrayList<Integer> removeHypothesis = new ArrayList<Integer>();

		    for (int objNum = 0; objNum < rsoChoices.size(); objNum++)
		    {
			for (int hypNum = 0; hypNum < rsoChoices.get(objNum).size(); hypNum++)
			{
			    final double stepStart = t0.durationFrom(odCfg.propStart);
			    final double stepEnd = tm.durationFrom(odCfg.propStart);
			    OneObject currSC = rsoChoices.get(objNum).get(hypNum);

			    if (currSC.dataAssociated)
			    {
				currSC.dataAssociated = false;
				currSC.sigma = generateSigmaPoints(currSC.xhat, currSC.P, numStates, numSigmas);
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
				catch(RuntimeException e)
				{
				    removeHypothesis.add(hypNum);
				    continue;
				}
			    }
			    else
				currSC.propSigma.setSubMatrix(currSC.sigma.getData(), 0, 0);

			    currSC.xhatPrev = addColumns(currSC.propSigma).mapMultiplyToSelf(weight);
			    currSC.xhat = new ArrayRealVector(currSC.xhatPrev);
			    currSC.Pprop = odCfg.getProcessNoiseMatrix(stepEnd - stepStart);
			    for (int i = 0; i < numSigmas; i++)
			    {
				RealVector y = currSC.propSigma.getColumnVector(i).subtract(currSC.xhatPrev);
				currSC.Pprop = currSC.Pprop.add(y.outerProduct(y).scalarMultiply(weight));
			    }

			    JPDALikelihoods JPDAtemp = new JPDALikelihoods();
			    JPDAtemp.psi = (1 - odCfg.estmDetectionProbability)*currSC.hypothesisWeight;
			    JPDAtemp.object = objNum;
			    JPDAtemp.measurement = 0;

			    currSC.marginalEvents.add(JPDAtemp);
			    // Compute predicted measurement, compare to each measurement					    
			    for (int measNum = 0; measNum <= additionalMeas; measNum++)
			    {
				RealVector rawMeas = updatePrep(currSC, tm, measIndex + measNum, numSigmas, propagator, biasPos);
				RealVector yhatpre = addColumns(currSC.estimMeas).mapMultiplyToSelf(weight);
				RealMatrix Pyy = R.copy();
				for (int i = 0; i < numSigmas; i++)
				{
				    RealVector y = currSC.estimMeas.getColumnVector(i).subtract(yhatpre);
				    Pyy = Pyy.add(y.outerProduct(y).scalarMultiply(weight));
				}

				RealVector Innov = rawMeas.subtract(yhatpre);
				RealMatrix MahalaTemp = MatrixUtils.createColumnRealMatrix(Innov.toArray());
				RealMatrix Mahalanobis = MahalaTemp.transpose().multiply(MatrixUtils.inverse(Pyy).multiply(MahalaTemp));
				if (Math.sqrt(Mahalanobis.getEntry(0, 0)) < odCfg.estmGatingThreshold)
				{
				    JPDAtemp = new JPDALikelihoods();
				    JPDAtemp.psi = (1 - new ChiSquaredDistribution(Rsize).cumulativeProbability(
							Mahalanobis.getEntry(0, 0)))*currSC.hypothesisWeight;
				    JPDAtemp.object = objNum;
				    JPDAtemp.measurement = measNum + 1;
				    currSC.marginalEvents.add(JPDAtemp);
				}
			    }
			}

			// Remove objects that fail to propagate
			for (int hypNum = 0; hypNum < removeHypothesis.size(); hypNum++)
			    rsoChoices.get(objNum).remove(removeHypothesis.get(hypNum));
		    }

		    // JPDA
		    double[][] sumJPDALikelihoods = new double[rsoChoices.size()][additionalMeas + 2]; // +1 to account for measurement 0 case
		    for (int objNum = 0; objNum < rsoChoices.size(); objNum++)
		    {
			for (int hypNum = 0; hypNum < rsoChoices.get(objNum).size(); hypNum++)
			{
			    OneObject currSC = rsoChoices.get(objNum).get(hypNum);
			    for (int eventCounter = 0; eventCounter < currSC.marginalEvents.size(); eventCounter++)
				sumJPDALikelihoods[objNum][currSC.marginalEvents.get(eventCounter).measurement] +=
				    currSC.marginalEvents.get(eventCounter).psi;
			}
		    }

		    jointEvents = JPDAJointEvents(jointEvents, sumJPDALikelihoods, singleJointEvent, 0);
		    double[] JPDAProbability = new double[jointEvents.size()];
		    Arrays.fill(JPDAProbability, 1.0);

		    for (int i = 0; i < jointEvents.size(); i++)
		    {
			for (int j = 0; j < jointEvents.get(i).size(); j++)
			{
			    double rowSum = 0;
			    for (int measNum = 0; measNum < sumJPDALikelihoods[j].length; measNum++)
				rowSum += sumJPDALikelihoods[j][measNum];
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
			    // Save Index
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
			    OneObject currSC = rsoChoices.get(objNum).get(hypNum);

			    // Update hypothesis weight
			    if (rsoChoices.get(objNum).size() > 1)
			    {
				double hypWeight = 0;			
				for (int counter = 0; counter < currSC.marginalEvents.size(); counter++)
				{
				    int measNum = currSC.marginalEvents.get(counter).measurement;
				    hypWeight += betaSatMeas[objNum][measNum]/sumJPDALikelihoods[objNum][measNum]*
					currSC.marginalEvents.get(counter).psi;
				}
				currSC.hypothesisWeight = hypWeight;
			    }

			    for (int i = 0; i < currSC.marginalEvents.size(); i++)
			    {
				if (currSC.marginalEvents.get(i).measurement > 0)
				{
				    double betaTemp = 0;
				    RealVector rawMeas = updatePrep(currSC, tm, measIndex + currSC.marginalEvents.get(i).measurement - 1,
								    numSigmas, propagator, biasPos);
				    RealVector yhatpre = addColumns(currSC.estimMeas).mapMultiplyToSelf(weight);
				    RealVector Innov = quadCheck(rawMeas, yhatpre);
				    if (currSC.marginalEvents.get(i).psi != 0)
				    {
					betaTemp = betaSatMeas[objNum][currSC.marginalEvents.get(i).measurement]/
					    sumJPDALikelihoods[objNum][currSC.marginalEvents.get(i).measurement]*
					    currSC.marginalEvents.get(i).psi/currSC.hypothesisWeight;
				    }
				    betav = betav.add(Innov.mapMultiply(betaTemp));
				    betavvt = betavvt.add(Innov.outerProduct(Innov).scalarMultiply(betaTemp));
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

			    currSC.xhat = new ArrayRealVector(currSC.xhatPrev.add(odCfg.parameterMatrix.multiply(K).operate(betav)));
			    currSC.P = currSC.Pprop.subtract(odCfg.parameterMatrix.multiply(K.multiply(Pyy.multiply(K.transpose()))).subtract(Ptilda));

			    // For Estimated Parameters, return xhat values to max/min bounds. If far outside the bounds, can result in
			    // Singular Matrix when taking the inverse for smoothing.
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
				Array2DRowRealMatrix postSigma = generateSigmaPoints(currSC.xhat, currSC.P, numStates, numSigmas);
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
				smout.measObjsSmoother = odObs.array[measIndex].helpers[0];
				currSC.smootherData.add(smout);
			    }
			}
		    }
		
		    // Merge/Prune hypotheses
		    for (int objNum = 0; objNum < rsoChoices.size(); objNum++)
		    {
			pruneHypotheses(rsoChoices.get(objNum));
			mergeHypotheses(rsoChoices.get(objNum));
		    }

		    // Normalize Weights
		    for (int objNum = 0; objNum < rsoChoices.size(); objNum++)
		    {
			double sumWeights = 0;
			for (int hypNum = 0; hypNum < rsoChoices.get(objNum).size(); hypNum++)
			    sumWeights += rsoChoices.get(objNum).get(hypNum).hypothesisWeight;
			for (int hypNum = 0; hypNum < rsoChoices.get(objNum).size(); hypNum++)
			    rsoChoices.get(objNum).get(hypNum).hypothesisWeight /= sumWeights;
		    }
		}

		// Smooth Data		
		for (int objNum = 0; objNum < rsoChoices.size(); objNum++)
		{
		    for (int hypNum = 0; hypNum < rsoChoices.get(objNum).size(); hypNum++)
		    {
			OneObject currSC = rsoChoices.get(objNum).get(hypNum);
			currSC.McReynoldsConsistencyPass = true;
			int smSize = currSC.smootherData.size() - 1;

			for (int i = 0; i < smSize; i++)
			{
			    SmootherTimeStep smDatak1 = currSC.smootherData.get(smSize - i);
			    SmootherTimeStep smDatak = currSC.smootherData.get(smSize - i - 1);
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

			// Compute McReynolds Consistency or merge all hypotheses into one hypothesis
			if (!rsoChoices.get(objNum).get(0).fromCAR)
			{
			    double McReynoldsVal = -99999;
			    for (int i = 0; i < currSC.smootherData.size()-1;i++)
			    {
				SmootherTimeStep smDatak1 = currSC.smootherData.get(smSize - i);
				SmootherTimeStep smDatak = currSC.smootherData.get(smSize - i - 1);
				RealMatrix delx = smDatak.xpost.subtract(smDatak.xstar);	
				RealMatrix delP = smDatak.Ppost.subtract(smDatak.Pstar);
				for (int j = 0; j < 5; j++)
				{			
				    McReynoldsVal = Math.max(McReynoldsVal, Math.abs(delx.getEntry(j,0))/Math.sqrt(Math.abs(delP.getEntry(j, j))));
				    if (Math.abs(delx.getEntry(j, 0)) / Math.sqrt(Math.abs(delP.getEntry(j, j))) >= 3)
					currSC.McReynoldsConsistencyPass = false;
				}
			    }

			    SpacecraftState[] smSsta = new SpacecraftState[1];		
			    for (int j = 0; j < currSC.smootherData.size();j++)
			    {
				double[] pv = currSC.smootherData.get(j).xstar.getColumn(0);
				tm = currSC.smootherData.get(j).tmSmoother;
				smSsta[0] = new SpacecraftState(new CartesianOrbit(new PVCoordinates(new Vector3D(pv[0], pv[1], pv[2]),
												     new Vector3D(pv[3], pv[4], pv[5])),
										   odCfg.propInertialFrame, tm, Constants.EGM96_EARTH_MU));

				// Compute smoothed residuals
				if (singleObject || measNames.length == 1)
				{
				    for (int i = 0; i < measNames.length; i++)
				    {		
					double[] fitv = currSC.smootherData.get(j).measObjsSmoother.estimate(0, 0, smSsta).getEstimatedValue();
					if (measNames.length == 1)
					    currSC.estOutput.get(j).postFit = fitv;
					else
					    currSC.estOutput.get(j).postFit = new double[]{fitv[i]};
				    }
				}
				else
				{
				    double[] fitv = currSC.smootherData.get(j).measObjsSmoother.estimate(0, 0, smSsta).getEstimatedValue();
				    currSC.estOutput.get(j).postFit = fitv;
				}

				currSC.estOutput.get(j).estimatedState = currSC.smootherData.get(j).xstar.getColumn(0);
				currSC.estOutput.get(j).estimatedCovariance = getLowerTriangle(currSC.smootherData.get(j).Pstar);
			    }
			}
			else
			{
			    rsoChoices.get(objNum).get(0).fromCAR = false;
			    currSC.McReynoldsConsistencyPass = false;
			    break;
			    // TODO Merge remaining hypotheses here
			}
		    }
		}

		// Check McReynolds and potentially break out
		boolean promotedTrackReset = false;
		for (Iterator<ArrayList<OneObject>> iter = rsoChoices.iterator(); iter.hasNext(); )
		{
		    OneObject currSC = iter.next().get(0);
		    if (currSC.McReynoldsConsistencyPass)
		    {
			promotedTrackReset = true;
			for (int idx = 0; idx < currSC.associatedObsIndex.size(); idx++)
			{
			    RealMatrix Ptemp = new Array2DRowRealMatrix(currSC.estOutput.get(idx).estimatedCovariance);
			    currSC.xhatPrev = new ArrayRealVector(currSC.estOutput.get(idx).estimatedState);
			    currSC.propSigma = generateSigmaPoints(currSC.xhatPrev, Ptemp, numStates, numSigmas);

			    RealMatrix Pyy = R.copy();
			    RealVector yhatpre = addColumns(currSC.estimMeas).mapMultiplyToSelf(weight);
			    RealVector rawMeas = updatePrep(currSC, tm, currSC.associatedObsIndex.get(idx), numSigmas, propagator, biasPos);
			    for (int i = 0; i < numSigmas; i++)
			    {
				RealVector y = currSC.estimMeas.getColumnVector(i).subtract(yhatpre);
				Pyy = Pyy.add(y.outerProduct(y).scalarMultiply(weight));
			    }
			    currSC.estOutput.get(idx).innovationCovariance = getLowerTriangle(Pyy);

			    // Add each measurement to the list of associated measurements.
			    removeObs.add(currSC.associatedObsIndex.get(idx));
			}

			promotedTracks.add(currSC);
			iter.remove();
		    }
		}

		if (promotedTrackReset)
		{
		    if (odObs.array.length - removeObs.size() <= 5)
			break;
		    if (rsoChoices.size() == 0)
		    {
			smoothIter = -1;
			activateCAR = true;
		    }
		}
	    }

	    // Save output
	    ArrayList<Integer> allAssoc = new ArrayList<Integer>(odObs.array.length);
	    multiOutput.estOutput = new ArrayList<ArrayList<Estimation.EstimationOutput>>(promotedTracks.size());
	    multiOutput.associatedObs = new ArrayList<ArrayList<Integer>>(odObs.array.length);
	    multiOutput.unassociatedObs = new ArrayList<Integer>(odObs.array.length);
	    for (OneObject obj: promotedTracks)
	    {
		multiOutput.estOutput.add(obj.estOutput);
		multiOutput.associatedObs.add(obj.associatedObsIndex);
		allAssoc.addAll(obj.associatedObsIndex);
	    }

	    for (int idx = 0; idx < odObs.array.length; idx++)
	    {
		if (!allAssoc.contains(idx))
		    multiOutput.unassociatedObs.add(idx);
	    }
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

	private void mergeHypotheses(ArrayList<OneObject> hypotheses)
	{
	    int hypNum1 = 0;
	    while (hypNum1 < hypotheses.size())
	    {	
		int hypNum2 = hypNum1 + 1;
		while (hypNum2 < hypotheses.size())
		{	
		    OneObject obj1 = hypotheses.get(hypNum1);
		    OneObject obj2 = hypotheses.get(hypNum2);
		    //Combine distributions
		    RealMatrix MahalaTemp = MatrixUtils.createColumnRealMatrix(obj1.xhat.subtract(obj2.xhat).toArray());
		    RealMatrix Mahalanobis1 = MahalaTemp.transpose().multiply(MatrixUtils.inverse(obj1.P).multiply(MahalaTemp));
		    RealMatrix Mahalanobis2 = MahalaTemp.transpose().multiply(MatrixUtils.inverse(obj2.P).multiply(MahalaTemp));

		    //If mahalanobis distance low, combine distributions.
		    if (Math.min(Mahalanobis1.getEntry(0, 0), Mahalanobis2.getEntry(0, 0)) < 1)
		    {
			RealVector xhatTemp =  obj1.xhat.mapMultiply(obj1.hypothesisWeight).add(obj2.xhat.mapMultiply(obj2.hypothesisWeight))
			    .mapMultiply(1/(obj1.hypothesisWeight + obj2.hypothesisWeight));
			RealMatrix matrixTemp = MatrixUtils.createColumnRealMatrix(obj1.xhat.toArray());
			RealMatrix PTemp = (obj1.P.add(matrixTemp.multiply(matrixTemp.transpose())))
			    .scalarMultiply(obj1.hypothesisWeight/(obj1.hypothesisWeight + obj2.hypothesisWeight));
			matrixTemp = MatrixUtils.createColumnRealMatrix(obj2.xhat.toArray());
			PTemp = PTemp.add((obj2.P.add(matrixTemp.multiply(matrixTemp.transpose())))
					  .scalarMultiply(obj2.hypothesisWeight/(obj1.hypothesisWeight + obj2.hypothesisWeight)));
			matrixTemp = MatrixUtils.createColumnRealMatrix(xhatTemp.toArray());
			hypotheses.get(hypNum1).xhat = new ArrayRealVector(xhatTemp);
			hypotheses.get(hypNum1).P = PTemp.subtract(matrixTemp.multiply(matrixTemp.transpose()));
			hypotheses.get(hypNum1).hypothesisWeight += hypotheses.get(hypNum2).hypothesisWeight;
			hypotheses.remove(hypNum2);
		    }
		    else
			hypNum2++;
		}
			
		hypNum1++;
	    }
	}

	private void pruneHypotheses(ArrayList<OneObject> hypotheses)
	{
	    double maxWeight = 0;
	    for (int hypNum = 0; hypNum < hypotheses.size(); hypNum++)
	    {
		if (hypotheses.get(hypNum).hypothesisWeight > maxWeight)
		    maxWeight = hypotheses.get(hypNum).hypothesisWeight;
	    }

	    int hypNum = 0;
	    while (hypNum < hypotheses.size())
	    {	
		if (hypotheses.get(hypNum).hypothesisWeight < maxWeight*1E-3)
		    hypotheses.remove(hypNum);
		else
		    hypNum++;
	    }
	}

	private	ArrayList<Hypothesis> generateOpticalHypotheses(Measurements.Measurement obs, double sigmaRA, double sigmaRAd, double sigmaDec,
								double sigmaDecd, int numStates, int numSigmas)
	{
	    AbsoluteDate date = obs.time;
	    double RA = obs.values[0];
	    double RA_d = obs.angleRates[0];
	    double Dec = obs.values[1];
	    double Dec_d = obs.angleRates[1];
	    TimeStampedPVCoordinates stationCoords = odCfg.stations.get(obs.station).getBaseFrame().getPVCoordinates(date, odCfg.propInertialFrame);
	    ArrayList <CAR.CARGaussianElement> CARGaussians = new CAR(RA, Dec, RA_d, Dec_d, stationCoords, 10000.0, 10.0, 10000.0 ,
								      17000000, 37000000, 0.05, singleObject).getCAR();
	    ArrayList<Hypothesis> objectHypotheses= new ArrayList<Hypothesis>();

	    for (int i = 0; i < CARGaussians.size(); i++)
	    {
		if (CARGaussians.get(i).ordinateStd < 0)  // Need to fix, These hypotheses should not occur in the first place.
		    continue;
		Hypothesis singleHypothesis= new Hypothesis();
		RealVector meanTemp = new ArrayRealVector(new double[] {CARGaussians.get(i).abscissaMean, RA, Dec,
			CARGaussians.get(i).ordinateMean, RA_d, Dec_d, odCfg.getInitialState()[6], odCfg.getInitialState()[7]});
		RealMatrix CovarTemp = new DiagonalMatrix(new double[] {Math.pow(CARGaussians.get(i).abscissaStd, 2), sigmaRA, sigmaDec,
			Math.pow(CARGaussians.get(i).ordinateStd, 2), sigmaRAd, sigmaDecd, Math.pow(odCfg.estmCovariance[6], 2),
			Math.pow(odCfg.estmCovariance[7], 2)});
		Array2DRowRealMatrix sigma = generateSigmaPoints(meanTemp, CovarTemp, numStates, numSigmas);
		Array2DRowRealMatrix sigmaXYZ = new Array2DRowRealMatrix(numStates, numSigmas);

		for (int j = 0; j < sigma.getColumnDimension(); j++)
		    sigma.setColumn(j, rangeRaDecToCart(sigma.getColumnVector(j), date, obs.station));
		new ManualPropagation(odCfg).propagate(0, sigma, CARGaussians.get(i).abscissaMean/Constants.SPEED_OF_LIGHT, sigmaXYZ, false);

		RealMatrix Pxx = new Array2DRowRealMatrix(numStates, numStates);
		RealVector xhat = addColumns(sigma).mapMultiplyToSelf(0.5/numStates);
		for (int j = 0; j < sigma.getColumnDimension(); j++)
		{
		    RealVector xhatdiff = sigma.getColumnVector(j).subtract(xhat);
		    Pxx = Pxx.add(xhatdiff.outerProduct(xhatdiff).scalarMultiply(0.5/numStates));
		}
		//return mean/covar structure
		singleHypothesis.xhat = new ArrayRealVector(xhat.toArray());
		singleHypothesis.P = Pxx;
		singleHypothesis.weight = CARGaussians.get(i).weight;
		objectHypotheses.add(singleHypothesis);
	    }
	    return(objectHypotheses);
	}

	private ArrayList<Hypothesis> generateRadarHypotheses(Measurements.Measurement obs, double sigmaRange, double sigmaRR,
							      double sigmaRA, double sigmaDec, int numStates, int numSigmas)
	{
	    AbsoluteDate date = obs.time;
	    double range = obs.values[0];
	    double rangeRate = obs.values[1];
	    double ra = obs.angleRates[0];
	    double dec = obs.angleRates[1];
	    TimeStampedPVCoordinates stationCoords = odCfg.stations.get(obs.station).getBaseFrame().getPVCoordinates(date, odCfg.propInertialFrame);

	    ArrayList <CAR.CARGaussianElement> CARGaussians = new CAR(ra, dec, range, rangeRate, stationCoords, 5e-6, 5e-6, 5e-6 ,
								      17000000, 37000000, 0.1, singleObject).getCAR();
	    ArrayList<Hypothesis> objectHypotheses= new ArrayList<Hypothesis>();

	    for (int i = 0; i < CARGaussians.size(); i++)
	    {
		if (CARGaussians.get(i).ordinateStd < 0)  //Need to fix, These hypotheses should not occur in the first place.
		    continue;
		Hypothesis singleHypothesis= new Hypothesis();
		RealVector meanTemp = new ArrayRealVector(new double[] {range, ra, dec, rangeRate, CARGaussians.get(i).abscissaMean,
			CARGaussians.get(i).ordinateMean, odCfg.getInitialState()[6], odCfg.getInitialState()[7]});
		RealMatrix CovarTemp = new DiagonalMatrix(new double[] {sigmaRange, sigmaRA, sigmaDec, sigmaRR,
			Math.pow(CARGaussians.get(i).abscissaStd, 2), Math.pow(CARGaussians.get(i).ordinateStd, 2),
			Math.pow(odCfg.estmCovariance[6], 2), Math.pow(odCfg.estmCovariance[7], 2)});
		Array2DRowRealMatrix sigma = generateSigmaPoints(meanTemp, CovarTemp, numStates, numSigmas);
		Array2DRowRealMatrix sigmaXYZ = new Array2DRowRealMatrix(numStates, numSigmas);

		for (int j = 0; j < sigma.getColumnDimension(); j++)
		    sigma.setColumn(j, rangeRaDecToCart(sigma.getColumnVector(j), date, obs.station));
		new ManualPropagation(odCfg).propagate(0, sigma, range/Constants.SPEED_OF_LIGHT, sigmaXYZ, false);

		RealMatrix Pxx = new Array2DRowRealMatrix(numStates, numStates);
		RealVector xhat = addColumns(sigma).mapMultiplyToSelf(0.5/numStates);
		for (int j = 0; j < sigma.getColumnDimension(); j++)
		{
		    RealVector xhatdiff = sigma.getColumnVector(j).subtract(xhat);
		    Pxx = Pxx.add(xhatdiff.outerProduct(xhatdiff).scalarMultiply(0.5/numStates));
		}
		//return mean/covar structure
		singleHypothesis.xhat = new ArrayRealVector(xhat.toArray());
		singleHypothesis.P = Pxx;
		singleHypothesis.weight = CARGaussians.get(i).weight;
		objectHypotheses.add(singleHypothesis);
	    }
	    return(objectHypotheses);
	}

	private	double[] rangeRaDecToCart(RealVector RangeRaDec, AbsoluteDate date, String stat)
	{
	    double Range = RangeRaDec.getEntry(0);
	    double Range_d = RangeRaDec.getEntry(3);		
	    double RA = RangeRaDec.getEntry(1);
	    double RA_d = RangeRaDec.getEntry(4);		
	    double Dec = RangeRaDec.getEntry(2);
	    double Dec_d = RangeRaDec.getEntry(5);

	    //Compute Inertial pos/vel relative to station
	    Vector3D topoPos = new Vector3D(new double[] {Range*Math.cos(Dec)*Math.cos(RA), Range*Math.cos(Dec)*Math.sin(RA), Range*Math.sin(Dec)});
	    Vector3D topoVel = new Vector3D(new double[] {Range_d*Math.cos(Dec)*Math.cos(RA) - Range*Math.sin(Dec)*Math.cos(RA)*Dec_d -
		    Range*Math.cos(Dec)*Math.sin(RA)*RA_d, Range_d*Math.cos(Dec)*Math.sin(RA) - Range*Math.sin(Dec)*Math.sin(RA)*Dec_d +
		    Range*Math.cos(Dec)*Math.cos(RA)*RA_d, Range_d*Math.sin(Dec) + Range*Math.cos(Dec) * Dec_d});

	    //get station coords
	    //add station coords to observation coords.
	    Vector3D Pos = new Vector3D(1, topoPos, 1, odCfg.stations.get(stat).getBaseFrame()
					.getPVCoordinates(date, odCfg.propInertialFrame).getPosition());
	    Vector3D Vel = new Vector3D(1, topoVel, 1, odCfg.stations.get(stat).getBaseFrame()
					.getPVCoordinates(date, odCfg.propInertialFrame).getVelocity());
	    return(new double[] {Pos.getX(), Pos.getY(), Pos.getZ(), Vel.getX(), Vel.getY(), Vel.getZ(),
		    RangeRaDec.getEntry(6), RangeRaDec.getEntry(7)});
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

	private RealVector quadCheck(RealVector measurement, RealVector state)
	{
	    //Innovations QuadCheck
	    RealVector Innov = measurement.subtract(state);
	    if (singleObject && Innov.getEntry(0) > Math.PI)
		Innov.setEntry(0, Innov.getEntry(0) - 2*Math.PI);
	    else if (singleObject && Innov.getEntry(0) < -Math.PI)
		Innov.setEntry(0, Innov.getEntry(0) + 2*Math.PI);
	    return(Innov);
	}

	private ArrayList<ArrayList<JPDALikelihoods>> JPDAJointEvents(ArrayList<ArrayList<JPDALikelihoods>> JointEvents, double[][] MarginalEvents,
								      ArrayList<JPDALikelihoods> SingleJointEvent, int objNum)
	{
	    for (int measNum = 0; measNum<MarginalEvents[objNum].length; measNum++)
	    {	
		ArrayList<JPDALikelihoods> SingleJointEventTemp = new ArrayList<JPDALikelihoods>();
		//Copy Single Joint Event into Single Joint Event Temp
		for (int m = 0; m < SingleJointEvent.size(); m++)
		    SingleJointEventTemp.add(SingleJointEvent.get(m));
		//Add event if that measurement has not been used
		boolean skip = false;
		for (int m = 0; m < SingleJointEvent.size(); m++)
		{
		    if (measNum == SingleJointEventTemp.get(m).measurement && measNum != 0)
			skip = true;
		}
	
		if (skip)
		    continue;
		JPDALikelihoods temp = new JPDALikelihoods();
		temp.measurement = measNum;
		temp.object = objNum;
		temp.psi = MarginalEvents[objNum][measNum];

		SingleJointEventTemp.add(temp);
		if (MarginalEvents.length == objNum + 1)
		    JointEvents.add(SingleJointEventTemp);
		else
		    JointEvents = JPDAJointEvents(JointEvents, MarginalEvents, SingleJointEventTemp, objNum+1);
	    }
	    return(JointEvents);
	}

	private Array2DRowRealMatrix generateSigmaPoints(RealVector mean, RealMatrix Cov, int numStates, int numSigmas)
	{
	    RealMatrix Ptemp;
	    if (Cov.getRowDimension() == Cov.getColumnDimension())
		Ptemp = Cov.scalarMultiply(numStates);
	    else
	    {
		Ptemp = new Array2DRowRealMatrix(numStates, numStates);
		for (int i = 0, k = 0; i < numStates; i++)
		{
		    for (int j = 0; j < i + 1; j++, k++)
		    {
			double value = Cov.getEntry(k, 0)*numStates;
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

	private RealVector updatePrep(OneObject currSC, AbsoluteDate tm, int measIndex, int numSigmas,
				      ManualPropagation propagator, HashMap<String, Integer> biasPos)
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

		if (singleObject || measNames.length == 1)
		{
		    double[] fitv = odObs.array[measIndex].helpers[0].estimate(0, 0, ssta).getEstimatedValue();
		    currSC.estimMeas.setColumn(i, fitv);
		    if (rawMeas == null)
			rawMeas = new ArrayRealVector(odObs.array[measIndex].helpers[0].getObservedValue());
		}
		else
		{
		    double[] fitv = odObs.array[measIndex].helpers[0].estimate(0, 0, ssta).getEstimatedValue();
		    currSC.estimMeas.setEntry(0, i, fitv[0]);
		    fitv = odObs.array[measIndex].helpers[1].estimate(0, 0, ssta).getEstimatedValue();
		    currSC.estimMeas.setEntry(1, i, fitv[0]);
		    if (rawMeas == null)
			rawMeas = new ArrayRealVector(new double[]{odObs.array[measIndex].helpers[0].getObservedValue()[0],
				odObs.array[measIndex].helpers[1].getObservedValue()[0]});
		}
	    }

	    //Does not enter if meas is pos/vel. Otherwise will enter.
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
    }
}
