/*
 * Estimation.java - Implementation of estimation algorithms.
 * Copyright (C) 2018-2022 University of Texas
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

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import org.hipparchus.geometry.euclidean.threed.Vector3D;
import org.hipparchus.linear.Array2DRowRealMatrix;
import org.hipparchus.linear.ArrayRealVector;
import org.hipparchus.linear.CholeskyDecomposition;
import org.hipparchus.linear.MatrixUtils;
import org.hipparchus.linear.RealMatrix;
import org.hipparchus.linear.RealVector;
import org.hipparchus.optim.nonlinear.vector.leastsquares.GaussNewtonOptimizer;
import org.hipparchus.util.FastMath;
import org.orekit.attitudes.AttitudeProvider;
import org.orekit.estimation.leastsquares.BatchLSEstimator;
import org.orekit.estimation.measurements.ObservableSatellite;
import org.orekit.estimation.measurements.Position;
import org.orekit.estimation.measurements.modifiers.OutlierFilter;
import org.orekit.estimation.sequential.CovarianceMatrixProvider;
import org.orekit.estimation.sequential.KalmanEstimation;
import org.orekit.estimation.sequential.KalmanEstimator;
import org.orekit.estimation.sequential.KalmanEstimatorBuilder;
import org.orekit.estimation.sequential.KalmanObserver;
import org.orekit.forces.ForceModel;
import org.orekit.orbits.CartesianOrbit;
import org.orekit.orbits.PositionAngle;
import org.orekit.propagation.Propagator;
import org.orekit.propagation.SpacecraftState;
import org.orekit.propagation.conversion.DormandPrince853IntegratorBuilder;
import org.orekit.propagation.conversion.NumericalPropagatorBuilder;
import org.orekit.propagation.sampling.OrekitFixedStepHandler;
import org.orekit.time.AbsoluteDate;
import org.orekit.utils.Constants;
import org.orekit.utils.ParameterDriversList;
import org.orekit.utils.PVCoordinates;

public final class Estimation
{
    public static class EstimationOutput
    {
        public AbsoluteDate time = AbsoluteDate.JULIAN_EPOCH;
        public String station = "";
        public double[] estimatedState;
        public double[] propagatedCovariance;
        public double[] innovationCovariance;
        public double[] estimatedCovariance;
        public double[] preFit;
        public double[] postFit;
        public Double clutterProbability;

        public EstimationOutput(AbsoluteDate time, String station)
        {
            this.time = time;
            this.station = station;
        }
    }

    private final Settings odCfg;
    private final Measurements odObs;

    private final int measSize;
    private final Settings.MeasurementType[] measNames;
    private final boolean singleObject;

    private ArrayList<EstimationOutput> estOutput;

    public final static String[] DMC_ACC_ESTM = {"DMC_ACC_X", "DMC_ACC_Y", "DMC_ACC_Z"};
    public final static String DMC_ACC_PROP = "DMC_STATE";

    public Estimation(Settings odCfg, Measurements odObs)
    {
        this.odCfg = odCfg;
        this.odObs = odObs;
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

    public ArrayList<EstimationOutput> determineOrbit()
    {
        estOutput = new ArrayList<EstimationOutput>(odObs.array.length);
        switch (odCfg.estmFilter)
        {
        case EXTENDED_KALMAN:
            new ExtendedKalmanFilter().determineOrbit();
            break;
        case UNSCENTED_KALMAN:
            new UnscentedKalmanFilter().determineOrbit();
            break;
        case BATCH_LEAST_SQUARES:
            new BatchLeastSquares().determineOrbit();
            break;
        default:
            throw(new RuntimeException("Invalid estimation filter"));
        }
        return(estOutput);
    }

    public static double[] getLowerTriangle(RealMatrix mat)
    {
        int m = mat.getRowDimension();
        double[] out = new double[(int)(0.5*m*(m+1))];
        for (int i = 0, k = 0; i < m; i++)
        {
            for (int j = 0; j <= i; j++)
                out[k++] = mat.getEntry(i, j);
        }
        return(out);
    }

    private class ExtendedKalmanFilter implements CovarianceMatrixProvider, KalmanObserver, OrekitFixedStepHandler
    {
        private int measIndex;
        private int stepIndex;
        private AbsoluteDate estmTime;
        private RealMatrix estmCovariance;

        private void determineOrbit()
        {
            final double[] x0 = odCfg.getInitialState();
            final CartesianOrbit orb0 = new CartesianOrbit(new PVCoordinates(new Vector3D(x0[0], x0[1], x0[2]), new Vector3D(x0[3], x0[4], x0[5])),
                                                           odCfg.propInertialFrame, odCfg.propStart, Constants.EGM96_EARTH_MU);
            estmTime = odCfg.propStart;
            estmCovariance = odCfg.getInitialCovariance();

            final PropagatorBuilder propBuilder = new PropagatorBuilder(odCfg, orb0, false);
            final KalmanEstimator filter = new KalmanEstimatorBuilder().addPropagationConfiguration(propBuilder, this).build();
            filter.setObserver(this);

            Propagator propagator = null;
            final Vector3D fakePosition = new Vector3D(Double.MAX_VALUE, Double.MAX_VALUE, Double.MAX_VALUE);
            final double[] fakeError = {100.0, 100.0, 100.0};
            final ObservableSatellite satellite = new ObservableSatellite(0);
            final OutlierFilter<Position> outlier = new OutlierFilter<Position>(1, 1.0);

            for (measIndex = 0; measIndex < odObs.array.length; measIndex++)
            {
                final Measurements.Measurement thisObs = odObs.array[measIndex];
                if (thisObs.values.length > 0)
                {
                    for (int i = 0; i < thisObs.helpers.length; i++)
                        propagator = filter.estimationStep(thisObs.helpers[i])[0];
                    propBuilder.enableDMC = true;
                }
                else if ((odCfg.outputFlags & Settings.OUTPUT_PROP_COV) != 0)
                {
                    Position fake = new Position(thisObs.time, fakePosition, fakeError, 1.0, satellite);
                    fake.addModifier(outlier);
                    filter.estimationStep(fake);
                }
                else
                {
                    int toIndex;
                    double stepSize;
                    if (measIndex > 0)
                        stepSize = thisObs.time.durationFrom(odObs.array[measIndex - 1].time);
                    else
                        stepSize = thisObs.time.durationFrom(odCfg.propStart);
                    for (toIndex = measIndex + 1; toIndex < odObs.array.length; toIndex++)
                    {
                        if (odObs.array[toIndex].values.length > 0 ||
                            odObs.array[toIndex].time.durationFrom(odObs.array[toIndex - 1].time) != stepSize)
                        {
                            toIndex--;
                            break;
                        }
                    }

                    stepIndex = measIndex;
                    if (propagator == null)
                        propagator = propBuilder.buildPropagator(propBuilder.getSelectedNormalizedParameters());
                    propagator.setStepHandler(stepSize, this);
                    propagator.propagate(odObs.array[toIndex].time);
                    measIndex = toIndex;
                }
            }
        }

        @Override public RealMatrix getInitialCovarianceMatrix(SpacecraftState init)
        {
            return(odCfg.getInitialCovariance());
        }

        @Override public RealMatrix getProcessNoiseMatrix(SpacecraftState prev, SpacecraftState curr)
        {
            return(odCfg.getProcessNoiseMatrix(curr.getDate().durationFrom(prev.getDate())));
        }

        @Override public void evaluationPerformed(KalmanEstimation est)
        {
            EstimationOutput result;
            final Measurements.Measurement thisObs = odObs.array[measIndex];
            if (measIndex >= estOutput.size())
            {
                result = new EstimationOutput(thisObs.time, thisObs.station);
                estOutput.add(result);
            }
            else
                result = estOutput.get(measIndex);

            final PVCoordinates pvc = est.getCorrectedSpacecraftStates()[0].getPVCoordinates();
            result.estimatedState = new double[odCfg.parameters.size() + 6];
            System.arraycopy(pvc.getPosition().toArray(), 0, result.estimatedState, 0, 3);
            System.arraycopy(pvc.getVelocity().toArray(), 0, result.estimatedState, 3, 3);

            ParameterDriversList plst = est.getEstimatedPropagationParameters();
            for (int i = 0; i < odCfg.parameters.size(); i++)
            {
                Settings.Parameter ep = odCfg.parameters.get(i);
                result.estimatedState[i + 6] = plst.findByName(ep.name).getValue();
            }

            if ((odCfg.outputFlags & Settings.OUTPUT_PROP_COV) != 0)
            {
                final RealMatrix phi = est.getPhysicalStateTransitionMatrix();
                if (estmCovariance != null && phi != null)
                    result.propagatedCovariance = getLowerTriangle(phi.multiply(estmCovariance).multiply(phi.transpose()).
                                                                   add(odCfg.getProcessNoiseMatrix(result.time.durationFrom(estmTime))));
                else if (est.getPhysicalEstimatedCovarianceMatrix() != null)
                    result.propagatedCovariance = getLowerTriangle(est.getPhysicalEstimatedCovarianceMatrix());
            }

            if (thisObs.values.length == 0)
                return;

            estmTime = result.time;
            estmCovariance = est.getPhysicalEstimatedCovarianceMatrix();
            if (estmCovariance != null && (odCfg.outputFlags & Settings.OUTPUT_ESTM_COV) != 0)
                result.estimatedCovariance = getLowerTriangle(estmCovariance);

            final double[] prev = est.getPredictedMeasurement().getEstimatedValue();
            final double[] post = est.getCorrectedMeasurement().getEstimatedValue();
            if (singleObject)
            {
                if ((odCfg.outputFlags & Settings.OUTPUT_RESIDUALS) != 0)
                {
                    result.preFit = prev;
                    result.postFit = post;
                }
                if (est.getPhysicalInnovationCovarianceMatrix() != null && (odCfg.outputFlags & Settings.OUTPUT_INNO_COV) != 0)
                    result.innovationCovariance = getLowerTriangle(est.getPhysicalInnovationCovarianceMatrix());
            }
            else
            {
                if ((odCfg.outputFlags & Settings.OUTPUT_RESIDUALS) != 0)
                {
                    if (result.preFit == null || result.postFit == null)
                    {
                        result.preFit = new double[measSize];
                        result.postFit = new double[measSize];
                        result.preFit[0] = prev[0];
                        result.postFit[0] = post[0];
                    }
                    else
                    {
                        result.preFit[1] = prev[0];
                        result.postFit[1] = post[0];
                    }
                }

                if (est.getPhysicalInnovationCovarianceMatrix() != null && (odCfg.outputFlags & Settings.OUTPUT_INNO_COV) != 0)
                {
                    if (result.innovationCovariance == null)
                    {
                        result.innovationCovariance = new double[(int)(0.5*measSize*(measSize+1))];
                        result.innovationCovariance[0] = est.getPhysicalInnovationCovarianceMatrix().getEntry(0, 0);
                    }
                    else
                        result.innovationCovariance[2] = est.getPhysicalInnovationCovarianceMatrix().getEntry(0, 0);
                }
            }
        }

        @Override public void handleStep(SpacecraftState state)
        {
            Measurements.Measurement thisObs = odObs.array[stepIndex];
            if (thisObs.time.equals(state.getDate()))
            {
                PVCoordinates pv = state.getPVCoordinates();
                EstimationOutput result = new EstimationOutput(thisObs.time, thisObs.station);
                if (stepIndex > 0)
                    result.estimatedState = estOutput.get(stepIndex - 1).estimatedState.clone();
                else
                    result.estimatedState = odCfg.getInitialState();
                System.arraycopy(pv.getPosition().toArray(), 0, result.estimatedState, 0, 3);
                System.arraycopy(pv.getVelocity().toArray(), 0, result.estimatedState, 3, 3);
                estOutput.add(result);
                stepIndex++;
            }
        }
    }

    private class UnscentedKalmanFilter
    {
        private void determineOrbit()
        {
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
                            if (bias.equals(odCfg.parameters.get(k).name))
                                biasPos.put(bias, k + 6);
                        }
                    }
                }
            }

            Array2DRowRealMatrix R = new Array2DRowRealMatrix(measSize, measSize);
            for (int i = 0, j = 0; i < measNames.length; i++)
            {
                Settings.Measurement jm = odCfg.cfgMeasurements.get(measNames[i]);
                for (int k = 0; k < jm.error.length; k++)
                {
                    R.setEntry(j, j, jm.error[k]*jm.error[k]);
                    j++;
                }
            }

            boolean enableDMC = false;
            final double[] xInitial = odCfg.getInitialState();
            AbsoluteDate startTime = odCfg.propStart;
            final int numStates = odCfg.parameters.size() + 6;
            final int numSigmas = 2*numStates;
            final double weight = 0.5/numStates;
            final SpacecraftState[] ssta = new SpacecraftState[1];
            final ManualPropagation propagator = new ManualPropagation(odCfg);
            RealMatrix P = odCfg.getInitialCovariance();
            final Array2DRowRealMatrix sigma = new Array2DRowRealMatrix(numStates, numSigmas);
            final Array2DRowRealMatrix propSigma = new Array2DRowRealMatrix(numStates, numSigmas);
            final Array2DRowRealMatrix estimMeas = new Array2DRowRealMatrix(measSize, numSigmas);
            final RealMatrix psdCorr = MatrixUtils.createRealIdentityMatrix(P.getRowDimension()).scalarMultiply(1.0E-6);
            RealMatrix Pprop = null;
            ArrayRealVector xhat = new ArrayRealVector(xInitial);
            RealVector xhatPrev = new ArrayRealVector(xInitial);
            final double[] noiseMult = new double[measSize];
            Arrays.fill(noiseMult, 1.0);

            for (int measIndex = 0; measIndex < odObs.array.length; measIndex++)
            {
                if (measIndex > 0)
                    startTime = odObs.array[measIndex - 1].time;
                final Measurements.Measurement thisObs = odObs.array[measIndex];
                final double stepStart = startTime.durationFrom(odCfg.propStart);
                final double stepFinal = thisObs.time.durationFrom(odCfg.propStart);
                final RealMatrix Ptemp = P.scalarMultiply(numStates);
                final RealMatrix sqrtP = new CholeskyDecomposition(
                    Ptemp.add(Ptemp.transpose()).scalarMultiply(0.5).add(psdCorr), 1E-6, 1E-16).getL();

                for (int i = 0; i < numStates; i++)
                {
                    sigma.setColumnVector(i, xhat.add(sqrtP.getColumnVector(i)));
                    sigma.setColumnVector(numStates + i, xhat.subtract(sqrtP.getColumnVector(i)));
                }

                double[][] sigData = sigma.getData();
                for (int j = 6; j < odCfg.parameters.size() + 6; j++)
                {
                    Settings.Parameter tempep = odCfg.parameters.get(j - 6);
                    for (int i = 0; i < numSigmas; i++)
                        sigData[j][i] = FastMath.min(FastMath.max(sigData[j][i], tempep.min), tempep.max);
                }
                sigma.setSubMatrix(sigData, 0, 0);

                if (stepFinal - stepStart > 1.0E-6)
                {
                    propagator.propagate(stepStart, sigma, stepFinal, propSigma, enableDMC);
                    sigma.setSubMatrix(propSigma.getData(), 0, 0);
                }
                else
                    propSigma.setSubMatrix(sigma.getData(), 0, 0);

                xhatPrev = addColumns(propSigma).mapMultiplyToSelf(weight);
                xhat = new ArrayRealVector(xhatPrev);
                Pprop = odCfg.getProcessNoiseMatrix(stepFinal - stepStart);
                for (int i = 0; i < numSigmas; i++)
                {
                    RealVector y = propSigma.getColumnVector(i).subtract(xhatPrev);
                    Pprop = Pprop.add(y.outerProduct(y).scalarMultiply(weight));
                }

                if (thisObs.values.length == 0)
                {
                    EstimationOutput odout = new EstimationOutput(thisObs.time, thisObs.station);
                    odout.estimatedState = xhatPrev.toArray();
                    if ((odCfg.outputFlags & Settings.OUTPUT_PROP_COV) != 0)
                        odout.propagatedCovariance = getLowerTriangle(Pprop);
                    estOutput.add(odout);
                    continue;
                }

                enableDMC = true;
                RealVector rawMeas = null, biasCorrection = null;
                for (int i = 0; i < numSigmas; i++)
                {
                    double[] pv = propSigma.getColumn(i);
                    CartesianOrbit orb = new CartesianOrbit(new PVCoordinates(new Vector3D(pv[0], pv[1], pv[2]), new Vector3D(pv[3], pv[4], pv[5])),
                                                            odCfg.propInertialFrame, thisObs.time, Constants.EGM96_EARTH_MU);
                    ssta[0] = new SpacecraftState(orb, odCfg.rsoMass);

                    if (singleObject)
                    {
                        double[] fitv = thisObs.helpers[0].estimate(0, 0, ssta).getEstimatedValue();
                        estimMeas.setColumn(i, fitv);
                        if (rawMeas == null)
                            rawMeas = new ArrayRealVector(thisObs.values);
                    }
                    else
                    {
                        double[] fitv = thisObs.helpers[0].estimate(0, 0, ssta).getEstimatedValue();
                        estimMeas.setEntry(0, i, fitv[0]);
                        fitv = thisObs.helpers[1].estimate(0, 0, ssta).getEstimatedValue();
                        estimMeas.setEntry(1, i, fitv[0]);
                        if (rawMeas == null)
                            rawMeas = new ArrayRealVector(thisObs.values);
                    }
                }

                if (thisObs.station.length() > 0)
                {
                    Integer pos = biasPos.get(new StringBuilder(thisObs.station).append(measNames[0]).toString());
                    if (pos != null)
                        biasCorrection = xhatPrev.getSubVector(pos, measSize);
                }

                for (int attempt = 0; attempt < 2; attempt++)
                {
                    RealVector yhatpre = addColumns(estimMeas).mapMultiplyToSelf(weight);
                    if (biasCorrection != null)
                        yhatpre = yhatpre.add(biasCorrection);

                    RealMatrix Pyy = R.copy();
                    for (int i = 0; attempt == 1 && i < measSize; i++)
                        Pyy.setEntry(i, i, Pyy.getEntry(i, i)*noiseMult[i]);

                    RealMatrix Pxy = new Array2DRowRealMatrix(numStates, measSize);
                    for (int i = 0; i < numSigmas; i++)
                    {
                        RealVector y = estimMeas.getColumnVector(i).subtract(yhatpre);
                        Pyy = Pyy.add(y.outerProduct(y).scalarMultiply(weight));
                        Pxy = Pxy.add(propSigma.getColumnVector(i).subtract(xhatPrev).outerProduct(y).scalarMultiply(weight));
                    }

                    EstimationOutput odout = new EstimationOutput(thisObs.time, thisObs.station);
                    RealVector error = rawMeas.subtract(yhatpre);
                    RealMatrix invPyy = MatrixUtils.inverse(Pyy);
                    RealMatrix K = Pxy.multiply(invPyy);
                    xhat = new ArrayRealVector(xhatPrev.add(odCfg.parameterMatrix.multiply(K).operate(error)));
                    P = Pprop.subtract(odCfg.parameterMatrix.multiply(K.multiply(Pyy.multiply(K.transpose()))));

                    double[] pv = xhat.toArray();
                    CartesianOrbit orb = new CartesianOrbit(new PVCoordinates(new Vector3D(pv[0], pv[1], pv[2]), new Vector3D(pv[3], pv[4], pv[5])),
                                                            odCfg.propInertialFrame, thisObs.time, Constants.EGM96_EARTH_MU);
                    ssta[0] = new SpacecraftState(orb, odCfg.rsoMass);

                    if ((odCfg.outputFlags & Settings.OUTPUT_RESIDUALS) != 0)
                    {
                        odout.preFit = yhatpre.toArray();
                        if (singleObject)
                            odout.postFit = thisObs.helpers[0].estimate(0, 0, ssta).getEstimatedValue();
                        else
                            odout.postFit = new double[]{thisObs.helpers[0].estimate(0, 0, ssta).getEstimatedValue()[0],
                                thisObs.helpers[1].estimate(0, 0, ssta).getEstimatedValue()[0]};
                    }

                    if (odCfg.estmOutlierSigma > 0.0 &&	odCfg.estmOutlierWarmup > 0 && measIndex >= odCfg.estmOutlierWarmup &&
                        attempt == 0 && odout.postFit != null)
                    {
                        int pos = 0;
                        boolean isOutlier = false;
                        double[] fit = odout.postFit;
                        for (int j = 0; j < fit.length; j++)
                        {
                            if (FastMath.pow(rawMeas.getEntry(pos) - fit[j], 2) >
                                odCfg.estmOutlierSigma*odCfg.estmOutlierSigma*Pyy.getEntry(pos, pos))
                            {
                                isOutlier = true;
                                noiseMult[pos] = 1.0E16;
                            }
                            else
                                noiseMult[pos] = 1.0;
                            pos++;
                        }

                        if (isOutlier)
                            continue;
                    }

                    estOutput.add(odout);
                    odout.estimatedState = pv;
                    if ((odCfg.outputFlags & Settings.OUTPUT_ESTM_COV) != 0)
                        odout.estimatedCovariance = getLowerTriangle(P);
                    if ((odCfg.outputFlags & Settings.OUTPUT_PROP_COV) != 0)
                        odout.propagatedCovariance = getLowerTriangle(Pprop);
                    if ((odCfg.outputFlags & Settings.OUTPUT_INNO_COV) != 0)
                        odout.innovationCovariance = getLowerTriangle(Pyy);
                    break;
                }
            }
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
    }

    private class BatchLeastSquares
    {
        private void determineOrbit()
        {
            double[] x0 = odCfg.getInitialState();
            CartesianOrbit X0 = new CartesianOrbit(new PVCoordinates(new Vector3D(x0[0], x0[1], x0[2]), new Vector3D(x0[3], x0[4], x0[5])),
                                                   odCfg.propInertialFrame, odCfg.propStart, Constants.EGM96_EARTH_MU);

            NumericalPropagatorBuilder propBuilder = new NumericalPropagatorBuilder(
                X0, new DormandPrince853IntegratorBuilder(odCfg.integMinTimeStep,odCfg.integMaxTimeStep, 1.0), PositionAngle.TRUE, 10.0);
            propBuilder.setMass(odCfg.rsoMass);
            for (ForceModel fm: odCfg.forces)
                propBuilder.addForceModel(fm);

            AttitudeProvider attProv = odCfg.getAttitudeProvider();
            if (attProv != null)
                propBuilder.setAttitudeProvider(attProv);

            BatchLSEstimator filter = new BatchLSEstimator(new GaussNewtonOptimizer(), propBuilder);
            filter.setParametersConvergenceThreshold(1E-3);
            filter.setMaxIterations(100);
            filter.setMaxEvaluations(300);

            for (Measurements.Measurement m: odObs.array)
            {
                for (int i = 0; i < m.helpers.length; i++)
                    filter.addMeasurement(m.helpers[i]);
            }

            Propagator propagator = filter.estimate()[0];
            for (Measurements.Measurement m: odObs.array)
            {
                EstimationOutput result = new EstimationOutput(m.time, m.station);
                estOutput.add(result);
                result.estimatedState = new double[6];
                PVCoordinates pv = propagator.getPVCoordinates(m.time, odCfg.propInertialFrame);
                System.arraycopy(pv.getPosition().toArray(), 0, result.estimatedState, 0, 3);
                System.arraycopy(pv.getVelocity().toArray(), 0, result.estimatedState, 3, 3);
            }
        }
    }
}
