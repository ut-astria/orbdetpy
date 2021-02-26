/*
 * Estimation.java - Implementation of estimation algorithms.
 * Copyright (C) 2018-2021 University of Texas
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
import java.util.List;
import org.hipparchus.geometry.euclidean.threed.Rotation;
import org.hipparchus.geometry.euclidean.threed.Vector3D;
import org.hipparchus.linear.Array2DRowRealMatrix;
import org.hipparchus.linear.ArrayRealVector;
import org.hipparchus.linear.CholeskyDecomposition;
import org.hipparchus.linear.DiagonalMatrix;
import org.hipparchus.linear.MatrixUtils;
import org.hipparchus.linear.RealMatrix;
import org.hipparchus.linear.RealVector;
import org.hipparchus.ode.ODEStateAndDerivative;
import org.hipparchus.util.FastMath;
import org.orekit.attitudes.AttitudeProvider;
import org.orekit.estimation.measurements.ObservedMeasurement;
import org.orekit.estimation.measurements.Range;
import org.orekit.estimation.sequential.ConstantProcessNoise;
import org.orekit.estimation.sequential.CovarianceMatrixProvider;
import org.orekit.estimation.sequential.KalmanEstimation;
import org.orekit.estimation.sequential.KalmanEstimator;
import org.orekit.estimation.sequential.KalmanEstimatorBuilder;
import org.orekit.estimation.sequential.KalmanObserver;
import org.orekit.forces.ForceModel;
import org.orekit.forces.gravity.NewtonianAttraction;
import org.orekit.orbits.CartesianOrbit;
import org.orekit.orbits.PositionAngle;
import org.orekit.propagation.SpacecraftState;
import org.orekit.propagation.conversion.DormandPrince853IntegratorBuilder;
import org.orekit.propagation.integration.AbstractIntegratedPropagator;
import org.orekit.propagation.numerical.NumericalPropagator;
import org.orekit.propagation.sampling.OrekitFixedStepHandler;
import org.orekit.time.AbsoluteDate;
import org.orekit.utils.Constants;
import org.orekit.utils.ParameterDriver;
import org.orekit.utils.ParameterDriversList;
import org.orekit.utils.PVCoordinates;
import org.orekit.utils.TimeStampedPVCoordinates;

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
    }

    private final Settings odCfg;
    private final Measurements odObs;

    private Settings.MeasurementType[] measNames;
    private final int measSize;

    private final AbsoluteDate epoch;
    private final AbsoluteDate propEnd;
    private final AbsoluteDate stepHandlerStart;
    private final AbsoluteDate stepHandlerEnd;
    private ArrayList<EstimationOutput> estOutput;

    public final static String DMC_ACC_ESTM[] = {"zDMCx", "zDMCy", "zDMCz"};
    public final static String DMC_ACC_PROP = "DMCEstProp";

    public Estimation(Settings odCfg, Measurements odObs)
    {
	this.odCfg = odCfg;
	this.odObs = odObs;
	if (odCfg.estmFilter == Settings.Filter.UNSCENTED_KALMAN && odCfg.gravityDegree >= 2 && odCfg.gravityOrder >= 0)
	    odCfg.forces.add(0, new NewtonianAttraction(Constants.EGM96_EARTH_MU));

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

	epoch = odCfg.propStart;
	propEnd = odCfg.propEnd;
	if (odCfg.propStep >= 0.0)
	{
	    stepHandlerStart = epoch;
	    stepHandlerEnd = propEnd;
	}
	else
	{
	    stepHandlerStart = propEnd;
	    stepHandlerEnd = epoch;
	}
    }

    public ArrayList<EstimationOutput> determineOrbit()
    {
	int size = FastMath.max(odObs.rawMeas.length, 10);
	if (odCfg.propStep != 0.0)
	    size += (int)FastMath.abs(propEnd.durationFrom(epoch)/odCfg.propStep) + 2;
	estOutput = new ArrayList<EstimationOutput>(size);

	if (odCfg.estmFilter == Settings.Filter.UNSCENTED_KALMAN)
	    new UnscentedKalmanFilter().determineOrbit();
	else
	    new ExtendedKalmanFilter().determineOrbit();
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
	private PropagatorBuilder propBuilder;
	private AbsoluteDate prevDate;
	private Vector3D prevPosition;
	private Vector3D prevVelocity;
	private RealMatrix prevCovariance;
	private double measDeltaT;
	private int measIndex;

	private void determineOrbit()
	{
	    final double[] Xi = odCfg.getInitialState();
	    prevDate = epoch;
	    prevPosition = new Vector3D(Xi[0],Xi[1],Xi[2]);
	    prevVelocity = new Vector3D(Xi[3],Xi[4],Xi[5]);
	    prevCovariance = odCfg.getInitialCovariance();
	    final CartesianOrbit X0 = new CartesianOrbit(new PVCoordinates(prevPosition, prevVelocity),
							 odCfg.propInertialFrame, epoch, Constants.EGM96_EARTH_MU);

	    OrekitFixedStepHandler handler = null;
    	    if (odCfg.propStep != 0.0)
		handler = this;

	    propBuilder = new PropagatorBuilder(odCfg, X0, new DormandPrince853IntegratorBuilder(odCfg.integMinTimeStep, odCfg.integMaxTimeStep, 1.0),
						PositionAngle.TRUE, 10.0, handler, false);
	    propBuilder.setMass(odCfg.rsoMass);
	    for (ForceModel fm: odCfg.forces)
		propBuilder.addForceModel(fm);

	    final ParameterDriversList plst = propBuilder.getPropagationParametersDrivers();
	    for (Settings.Parameter ep : odCfg.parameters)
	    {
		ParameterDriver pdrv = new ParameterDriver(ep.name, ep.value, 1.0, ep.min, ep.max);
		pdrv.setSelected(true);
		plst.add(pdrv);
	    }

	    final AttitudeProvider attprov = odCfg.getAttitudeProvider();
	    if (attprov != null)
		propBuilder.setAttitudeProvider(attprov);

	    final KalmanEstimatorBuilder builder = new KalmanEstimatorBuilder();
	    builder.addPropagationConfiguration(propBuilder, this);
	    final KalmanEstimator filter = builder.build();
	    filter.setObserver(this);

	    AbstractIntegratedPropagator estimator = propBuilder.buildPropagator(propBuilder.getSelectedNormalizedParameters());
	    for (measIndex = 0; measIndex < odObs.rawMeas.length; measIndex++)
	    {
		ObservedMeasurement m0 = null, m1 = null;
		if (measNames[0] != Settings.MeasurementType.RANGE && measNames[0] != Settings.MeasurementType.RANGE_RATE || measNames.length == 1)
		{
		    m0 = odObs.measObjs.get(measIndex);
		}
		else
		{
		    m0 = odObs.measObjs.get(measIndex*2);
		    m1 = odObs.measObjs.get(measIndex*2 + 1);
		}

		if (m0.isEnabled())
		{
		    AbstractIntegratedPropagator[] estimators = filter.estimationStep(m0);
		    if (m1 != null && m1.isEnabled())
			estimators = filter.estimationStep(m1);
		    if (estimators != null)
			estimator = estimators[0];
		}
		else
		    handleStep(estimator.propagate(m0.getDate()), false);
	    }

	    propBuilder.enableDMC = false;
	    if (handler == null && (odObs.rawMeas.length == 0 || !odCfg.propEnd.equals(odObs.rawMeas[odObs.rawMeas.length - 1].time)))
		handleStep(estimator.propagate(propEnd), true);
	}

	@Override public RealMatrix getInitialCovarianceMatrix(SpacecraftState init)
	{
	    return(odCfg.getInitialCovariance());
	}

	@Override public RealMatrix getProcessNoiseMatrix(SpacecraftState prev, SpacecraftState curr)
	{
	    double tmeas = curr.getDate().durationFrom(prev.getDate());
	    if (tmeas != 0.0)
		measDeltaT = tmeas;
	    else
		tmeas = measDeltaT;
	    return(odCfg.getProcessNoiseMatrix(tmeas));
	}

	@Override public void evaluationPerformed(KalmanEstimation est)
	{
	    propBuilder.enableDMC = true;
	    EstimationOutput result = null;
	    for (EstimationOutput loop: estOutput)
	    {
		if (loop.time.equals(odObs.rawMeas[measIndex].time) && loop.station.equals(odObs.rawMeas[measIndex].station))
		{
		    result = loop;
		    break;
		}
	    }

	    if (result == null)
	    {
		result = new EstimationOutput();
		result.time = odObs.rawMeas[measIndex].time;
		result.station = odObs.rawMeas[measIndex].station;
		estOutput.add(result);
	    }

	    final SpacecraftState ssta = est.getCorrectedSpacecraftStates()[0];
	    final PVCoordinates pvc = ssta.getPVCoordinates();
	    result.estimatedState = new double[odCfg.parameters.size() + 6];
	    System.arraycopy(pvc.getPosition().toArray(), 0, result.estimatedState, 0, 3);
	    System.arraycopy(pvc.getVelocity().toArray(), 0, result.estimatedState, 3, 3);

	    final ParameterDriversList plst = est.getEstimatedPropagationParameters();
	    for (int i = 0; i < odCfg.parameters.size(); i++)
	    {
		Settings.Parameter ep = odCfg.parameters.get(i);
		result.estimatedState[i + 6] = est.getEstimatedPropagationParameters().findByName(ep.name).getValue();
	    }

	    if (est.getPhysicalEstimatedCovarianceMatrix() != null && (odCfg.outputFlags & Settings.OUTPUT_ESTM_COV) != 0)
		result.estimatedCovariance = getLowerTriangle(est.getPhysicalEstimatedCovarianceMatrix());

	    final double[] prev = est.getPredictedMeasurement().getEstimatedValue();
	    final double[] post = est.getCorrectedMeasurement().getEstimatedValue();
	    if (measNames[0] != Settings.MeasurementType.RANGE && measNames[0] != Settings.MeasurementType.RANGE_RATE)
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

	    RealMatrix phi = est.getPhysicalStateTransitionMatrix();
	    if (prevCovariance != null && phi != null && (odCfg.outputFlags & Settings.OUTPUT_PROP_COV) != 0)
	    {
		result.propagatedCovariance = getLowerTriangle(
		    phi.multiply(prevCovariance).multiply(phi.transpose()).add(odCfg.getProcessNoiseMatrix(measDeltaT)));
		prevCovariance = est.getPhysicalEstimatedCovarianceMatrix();
	    }
	}

	@Override public void handleStep(SpacecraftState state, boolean lastStep)
	{
	    final PVCoordinates pvc = state.getPVCoordinates(odCfg.propInertialFrame);
	    for (ObservedMeasurement m: odObs.measObjs)
	    {
		if (m.isEnabled() && m.getDate().durationFrom(state.getDate()) == 0.0)
		    return;
	    }

	    final EstimationOutput odout = new EstimationOutput();
	    estOutput.add(odout);
	    odout.time = state.getDate();
	    odout.estimatedState = new double[odCfg.parameters.size() + 6];
	    System.arraycopy(pvc.getPosition().toArray(), 0, odout.estimatedState, 0, 3);
	    System.arraycopy(pvc.getVelocity().toArray(), 0, odout.estimatedState, 3, 3);
	    if (odCfg.parameters.size() > 0 && estOutput.size() > 0)
		System.arraycopy(estOutput.get(estOutput.size()-1).estimatedState, 6, odout.estimatedState, 6, odCfg.parameters.size());

	    if ((odCfg.outputFlags & Settings.OUTPUT_PROP_COV) != 0 && prevCovariance != null)
	    {
		Rotation phi1 = new Rotation(prevPosition, pvc.getPosition());
		Rotation phi2 = new Rotation(prevVelocity, pvc.getVelocity());
		RealMatrix phi = MatrixUtils.createRealIdentityMatrix(prevCovariance.getRowDimension());
		phi.setSubMatrix(phi1.getMatrix(), 0, 0);
		phi.setSubMatrix(phi2.getMatrix(), 3, 3);
		prevPosition = pvc.getPosition();
		prevVelocity = pvc.getVelocity();
		odout.propagatedCovariance = getLowerTriangle(phi.multiply(prevCovariance).multiply(phi.transpose()).add(
								  odCfg.getProcessNoiseMatrix(state.getDate().durationFrom(prevDate))));
		if (odObs.rawMeas.length > 0)
		    prevDate = state.getDate();
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
	    final int numStates = odCfg.parameters.size() + 6;
	    final int numSigmas = 2*numStates;
	    final double weight = 0.5/numStates;
	    final double[] xInitial = odCfg.getInitialState();
	    final double bound0 = stepHandlerStart.durationFrom(epoch);
	    final double bound1 = stepHandlerEnd.durationFrom(epoch);
	    AbsoluteDate tm = epoch;
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

	    for (int measIndex = 0; measIndex <= odObs.rawMeas.length; measIndex++)
	    {
		final AbsoluteDate t0 = tm;
		if (measIndex < odObs.rawMeas.length)
		    tm = odObs.rawMeas[measIndex].time;
		else
		{
		    tm = propEnd;
		    enableDMC = false;
		    if (FastMath.abs(tm.durationFrom(t0)) <= 1E-6)
			break;
		}

		double stepStart = t0.durationFrom(epoch), stepSum = 0.0;
		final double stepFinal = tm.durationFrom(epoch);
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

		boolean useObs = false;
		for (double v: odObs.rawMeas[measIndex].values)
		{
		    if (v != 0.0)
		    {
			useObs = true;
			break;
		    }
		}

		while (true)
		{
		    double step = stepFinal - stepStart;
		    if (odCfg.propStep != 0.0 && (stepStart >= bound0 || stepFinal >= bound0) &&
			(stepStart <= bound1 || stepFinal <= bound1))
		    {
			if (odCfg.propStep > 0.0)
			    step = FastMath.min(step, odCfg.propStep);
			else
			    step = FastMath.max(step, odCfg.propStep);
		    }
		    stepSum += step;

		    if (FastMath.abs(step) > 1.0E-6)
		    {
			propagator.propagate(stepStart, sigma, stepStart + step, propSigma, enableDMC);
			sigma.setSubMatrix(propSigma.getData(), 0, 0);
		    }
		    else
			propSigma.setSubMatrix(sigma.getData(), 0, 0);

		    xhatPrev = addColumns(propSigma).mapMultiplyToSelf(weight);
		    xhat = new ArrayRealVector(xhatPrev);
		    Pprop = odCfg.getProcessNoiseMatrix(stepSum);
		    for (int i = 0; i < numSigmas; i++)
		    {
			RealVector y = propSigma.getColumnVector(i).subtract(xhatPrev);
			Pprop = Pprop.add(y.outerProduct(y).scalarMultiply(weight));
		    }

		    if (measIndex == odObs.rawMeas.length || !useObs ||
			(odCfg.propStep != 0.0 && stepStart + step >= bound0 && stepStart + step <= bound1))
		    {
			EstimationOutput odout = new EstimationOutput();
			odout.time = new AbsoluteDate(epoch, stepStart + step);
			odout.estimatedState = xhatPrev.toArray();
			if ((odCfg.outputFlags & Settings.OUTPUT_PROP_COV) != 0)
			    odout.propagatedCovariance = getLowerTriangle(Pprop);
			estOutput.add(odout);
		    }

		    stepStart += step;
		    if (FastMath.abs(step) < 1.0E-6 || FastMath.abs(stepFinal - stepStart) < 1.0E-6)
			break;
		}

		if (measIndex == odObs.rawMeas.length)
		    break;
		if (!useObs)
		    continue;

		enableDMC = true;
		RealVector rawMeas = null, biasCorrection = null;
		for (int i = 0; i < numSigmas; i++)
		{
		    double[] pv = propSigma.getColumn(i);
		    ssta[0] = new SpacecraftState(new CartesianOrbit(new PVCoordinates(new Vector3D(pv[0], pv[1], pv[2]),
										       new Vector3D(pv[3], pv[4], pv[5])),
								     odCfg.propInertialFrame, tm, Constants.EGM96_EARTH_MU),
						  propagator.getAttitude(tm, pv), odCfg.rsoMass);

		    if (measNames[0] != Settings.MeasurementType.RANGE && measNames[0] != Settings.MeasurementType.RANGE_RATE ||
			measNames.length == 1)
		    {
			double[] fitv = odObs.measObjs.get(measIndex).estimate(0, 0, ssta).getEstimatedValue();
			estimMeas.setColumn(i, fitv);
			if (rawMeas == null)
			    rawMeas = new ArrayRealVector(odObs.measObjs.get(measIndex).getObservedValue());
		    }
		    else
		    {
			double[] fitv = odObs.measObjs.get(measIndex*2).estimate(0, 0, ssta).getEstimatedValue();
			estimMeas.setEntry(0, i, fitv[0]);
			fitv = odObs.measObjs.get(measIndex*2 + 1).estimate(0, 0, ssta).getEstimatedValue();
			estimMeas.setEntry(1, i, fitv[0]);
			if (rawMeas == null)
			    rawMeas = new ArrayRealVector(new double[]{odObs.measObjs.get(measIndex*2).getObservedValue()[0],
								       odObs.measObjs.get(measIndex*2 + 1).getObservedValue()[0]});
		    }
		}

		if (odObs.rawMeas[measIndex].station.length() > 0)
		{
		    Integer pos = biasPos.get(new StringBuilder(odObs.rawMeas[measIndex].station).append(measNames[0]).toString());
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

		    EstimationOutput odout = new EstimationOutput();
		    RealVector error = rawMeas.subtract(yhatpre);
		    RealMatrix invPyy = MatrixUtils.inverse(Pyy);
		    RealMatrix K = Pxy.multiply(invPyy);
		    P = Pprop.subtract(odCfg.parameterMatrix.multiply(K.multiply(Pyy.multiply(K.transpose()))));

		    if (odCfg.estmEnablePDAF)
		    {
			double alpha = FastMath.exp(-0.5*error.dotProduct(invPyy.operate(error)));
			double b = 2.0*(1.0-odCfg.estmGatingProbability*odCfg.estmDetectionProbability)*
			    FastMath.sqrt(new CholeskyDecomposition(Pyy, 1E-6, 1E-16).getDeterminant())/
			    (odCfg.estmDetectionProbability*odCfg.estmGatingThreshold);
			odout.clutterProbability = b/(b + alpha);
			xhat = new ArrayRealVector(xhatPrev.add(odCfg.parameterMatrix.multiply(K)
								.operate(error.mapMultiply(1-odout.clutterProbability))));
			RealMatrix yyt = error.outerProduct(error).scalarMultiply(
			    odout.clutterProbability*(1-odout.clutterProbability));
			P = Pprop.scalarMultiply(odout.clutterProbability).add(
			    P.scalarMultiply(1-odout.clutterProbability)).add(
				odCfg.parameterMatrix.multiply(K.multiply(yyt.multiply(K.transpose()))));
		    }
		    else
			xhat = new ArrayRealVector(xhatPrev.add(odCfg.parameterMatrix.multiply(K).operate(error)));

		    double[] pv = xhat.toArray();
		    ssta[0] = new SpacecraftState(new CartesianOrbit(new PVCoordinates(new Vector3D(pv[0], pv[1], pv[2]),
										       new Vector3D(pv[3], pv[4], pv[5])),
								     odCfg.propInertialFrame, tm, Constants.EGM96_EARTH_MU),
						  propagator.getAttitude(tm, pv), odCfg.rsoMass);

		    if ((odCfg.outputFlags & Settings.OUTPUT_RESIDUALS) != 0)
		    {
			odout.preFit = yhatpre.toArray();
			if (measNames[0] != Settings.MeasurementType.RANGE && measNames[0] != Settings.MeasurementType.RANGE_RATE ||
			    measNames.length == 1)
			    odout.postFit = odObs.measObjs.get(measIndex).estimate(0, 0, ssta).getEstimatedValue();
			else
			    odout.postFit = new double[]{odObs.measObjs.get(measIndex*2).estimate(0, 0, ssta).getEstimatedValue()[0],
							 odObs.measObjs.get(measIndex*2+1).estimate(0, 0, ssta).getEstimatedValue()[0]};
		    }

		    if (odCfg.estmOutlierSigma > 0.0 &&	odCfg.estmOutlierWarmup > 0 && measIndex >= odCfg.estmOutlierWarmup &&
			!odCfg.estmEnablePDAF && attempt == 0 && odout.postFit != null)
		    {
			int pos = 0;
			boolean isOutlier = false;
			for (int i = 0; i < measNames.length; i++)
			{
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
			}

			if (isOutlier)
			    continue;
		    }

		    estOutput.add(odout);
		    odout.time = odObs.rawMeas[measIndex].time;
		    odout.station = odObs.rawMeas[measIndex].station;
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
}
