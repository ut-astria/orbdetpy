/*
 * Estimation.java - Implementation of estimation algorithms.
 * Copyright (C) 2018-2019 University of Texas
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
	public String time;
	public String station;
	public double[] estimatedState;
	public double[][] propagatedCovariance;
	public double[][] innovationCovariance;
	public double[][] estimatedCovariance;
	public HashMap<String, double[]> preFit;
	public HashMap<String, double[]> postFit;
    }

    private final Settings odCfg;
    private final Measurements odObs;

    private String[] measNames;
    private final boolean combinedMeas;

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

	if (odCfg.estmFilter.equalsIgnoreCase("UKF") && odCfg.gravityDegree >= 2 && odCfg.gravityOrder >= 0)
	    odCfg.forces.add(0, new NewtonianAttraction(Constants.EGM96_EARTH_MU));

	measNames = odCfg.cfgMeasurements.keySet().toArray(new String[0]);
	Arrays.sort(measNames);
	if (measNames[0].equalsIgnoreCase("Declination"))
	    measNames = new String[]{"RightAscension", "Declination"};
	combinedMeas = !measNames[0].equalsIgnoreCase("Range") && !measNames[0].equalsIgnoreCase("RangeRate");

	epoch = DataManager.parseDateTime(odCfg.propStart);
	propEnd = DataManager.parseDateTime(odCfg.propEnd);

	if (odCfg.propStepHandlerStartTime != null)
	    stepHandlerStart = DataManager.parseDateTime(odCfg.propStepHandlerStartTime);
	else
	{
	    if (odCfg.propStep >= 0.0)
		stepHandlerStart = epoch;
	    else
		stepHandlerStart = propEnd;
	}

	if (odCfg.propStepHandlerEndTime != null)
	    stepHandlerEnd = DataManager.parseDateTime(odCfg.propStepHandlerEndTime);
	else
	{
	    if (odCfg.propStep >= 0.0)
		stepHandlerEnd = propEnd;
	    else
		stepHandlerEnd = epoch;
	}
    }

    public ArrayList<EstimationOutput> determineOrbit()
    {
	int size = FastMath.max(odObs.rawMeas.length, 10);
	if (odCfg.propStep != 0.0)
	    size += (int) FastMath.abs(propEnd.durationFrom(epoch)/odCfg.propStep) + 2;
	estOutput = new ArrayList<EstimationOutput>(size);

	if (odCfg.estmFilter.equalsIgnoreCase("UKF"))
	    new UnscentedKalmanFilter().determineOrbit();
	else
	    new ExtendedKalmanFilter().determineOrbit();

	return(estOutput);
    }

    private class ExtendedKalmanFilter implements CovarianceMatrixProvider, KalmanObserver, OrekitFixedStepHandler
    {
	private PropagatorBuilder propBuilder;
	private AbsoluteDate prevDate;
	private Vector3D prevPosition;
	private Vector3D prevVelocity;
	private RealMatrix prevCovariance;
	private double measDeltaT;

	private void determineOrbit()
	{
	    final double[] Xi = odCfg.getInitialState();
	    prevDate = epoch;
	    prevPosition = new Vector3D(Xi[0],Xi[1],Xi[2]);
	    prevVelocity = new Vector3D(Xi[3],Xi[4],Xi[5]);
	    prevCovariance = new DiagonalMatrix(odCfg.estmCovariance);
	    final CartesianOrbit X0 = new CartesianOrbit(new PVCoordinates(prevPosition, prevVelocity),
							 odCfg.propFrame, epoch, Constants.EGM96_EARTH_MU);

	    OrekitFixedStepHandler handler = null;
    	    if (odCfg.propStep != 0.0)
		handler = this;

	    propBuilder = new PropagatorBuilder(odCfg, X0, new DormandPrince853IntegratorBuilder(
						    odCfg.integMinTimeStep, odCfg.integMaxTimeStep, 1.0),
						PositionAngle.TRUE, 10.0, handler, false);
	    propBuilder.setMass(odCfg.rsoMass);
	    for (ForceModel fm : odCfg.forces)
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

	    AbstractIntegratedPropagator estimator = null;
	    AbstractIntegratedPropagator[] estimators = filter.processMeasurements(odObs.measObjs);
	    propBuilder.enableDMC = false;
	    if (estimators != null)
		estimator = estimators[0];
	    else
		estimator = propBuilder.buildPropagator(propBuilder.getSelectedNormalizedParameters());

	    if (odObs.rawMeas.length == 0 || !odCfg.propEnd.equalsIgnoreCase(odObs.rawMeas[odObs.rawMeas.length-1].time))
	    {
		final SpacecraftState state = estimator.propagate(propEnd);
		if (handler == null)
		    handleStep(state, true);
	    }
	}

	@Override public RealMatrix getInitialCovarianceMatrix(SpacecraftState init)
	{
	    return(new DiagonalMatrix(odCfg.estmCovariance));
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
	    int n = est.getCurrentMeasurementNumber() - 1;
	    if (!combinedMeas)
		n /= measNames.length;
	    final Measurements.Measurement raw = odObs.rawMeas[n];

	    EstimationOutput res = null;
	    for (EstimationOutput loop: estOutput)
	    {
		if (loop.preFit != null && loop.postFit != null && loop.time.equalsIgnoreCase(raw.time) &&
		    (loop.station == null || loop.station.equalsIgnoreCase(raw.station)))
		{
		    res = loop;
		    break;
		}
	    }

	    String key;
	    if (res == null)
	    {
		key = measNames[0];
		res = new EstimationOutput();
		res.preFit = new HashMap<String, double[]>();
		res.postFit = new HashMap<String, double[]>();
		res.time = odObs.rawMeas[n].time;
		res.station = odObs.rawMeas[n].station;
		estOutput.add(res);
	    }
	    else
		key = measNames[1];

	    final SpacecraftState ssta = est.getPredictedSpacecraftStates()[0];
	    final PVCoordinates pvc = ssta.getPVCoordinates();
	    res.estimatedState = new double[odCfg.parameters.size() + 6];
	    System.arraycopy(pvc.getPosition().toArray(), 0, res.estimatedState, 0, 3);
	    System.arraycopy(pvc.getVelocity().toArray(), 0, res.estimatedState, 3, 3);

	    final ParameterDriversList plst = est.getEstimatedPropagationParameters();
	    for (int i = 0; i < odCfg.parameters.size(); i++)
	    {
		Settings.Parameter ep = odCfg.parameters.get(i);
		res.estimatedState[i + 6] = est.getEstimatedPropagationParameters().findByName(ep.name).getValue();
	    }

	    final double[] pre = est.getPredictedMeasurement().getEstimatedValue();
	    final double[] pos = est.getCorrectedMeasurement().getEstimatedValue();
	    if (combinedMeas)
	    {
		for (int i = 0; i < measNames.length; i++)
		{
		    if (measNames.length == 1)
		    {
			res.preFit.put(measNames[i], pre);
			res.postFit.put(measNames[i], pos);
		    }
		    else
		    {
			res.preFit.put(measNames[i], new double[] {pre[i]});
			res.postFit.put(measNames[i], new double[] {pos[i]});
		    }
		}

		if (est.getPhysicalInnovationCovarianceMatrix() != null)
		    res.innovationCovariance = est.getPhysicalInnovationCovarianceMatrix().getData();
	    }
	    else
	    {
		res.preFit.put(key, pre);
		res.postFit.put(key, pos);
		if (est.getPhysicalInnovationCovarianceMatrix() != null)
		{
		    if (res.innovationCovariance == null)
			res.innovationCovariance = new double[][]{{est.getPhysicalInnovationCovarianceMatrix().getData()[0][0], 0.0}, {0.0, 0.0}};
		    else
			res.innovationCovariance[1][1] = est.getPhysicalInnovationCovarianceMatrix().getData()[0][0];
		}
	    }

	    final RealMatrix phi = est.getPhysicalStateTransitionMatrix();
	    if (phi != null && prevCovariance != null)
	    {
		res.propagatedCovariance = phi.multiply(prevCovariance).multiply(phi.transpose()).add(
		    odCfg.getProcessNoiseMatrix(measDeltaT)).getData();
	    }

	    if (est.getPhysicalEstimatedCovarianceMatrix() != null)
		res.estimatedCovariance = est.getPhysicalEstimatedCovarianceMatrix().getData();
	    prevCovariance = est.getPhysicalEstimatedCovarianceMatrix();
	}

	@Override public void handleStep(SpacecraftState state, boolean lastStep)
	{
	    if (state.getDate().durationFrom(stepHandlerStart) < 0.0 || state.getDate().durationFrom(stepHandlerEnd) > 0.0)
		return;

	    final PVCoordinates pvc = state.getPVCoordinates(odCfg.propFrame);
	    for (ObservedMeasurement m: odObs.measObjs)
	    {
		if (m.getDate().durationFrom(state.getDate()) == 0.0)
		    return;
	    }

	    final EstimationOutput odout = new EstimationOutput();
	    odout.time = DataManager.getUTCString(state.getDate());
	    odout.estimatedState = new double[odCfg.parameters.size() + 6];
	    System.arraycopy(pvc.getPosition().toArray(), 0, odout.estimatedState, 0, 3);
	    System.arraycopy(pvc.getVelocity().toArray(), 0, odout.estimatedState, 3, 3);
	    if (odCfg.parameters.size() > 0 && estOutput.size() > 0)
	    {
		System.arraycopy(estOutput.get(estOutput.size() - 1).estimatedState, 6,
				 odout.estimatedState, 6, odCfg.parameters.size());
	    }

	    final Rotation phi1 = new Rotation(prevPosition, pvc.getPosition());
	    final Rotation phi2 = new Rotation(prevVelocity, pvc.getVelocity());
	    final RealMatrix phi = MatrixUtils.createRealIdentityMatrix(prevCovariance.getRowDimension());
	    phi.setSubMatrix(phi1.getMatrix(), 0, 0);
	    phi.setSubMatrix(phi2.getMatrix(), 3, 3);
	    prevPosition = pvc.getPosition();
	    prevVelocity = pvc.getVelocity();

	    odout.propagatedCovariance = phi.multiply(prevCovariance).multiply(phi.transpose()).add(
		odCfg.getProcessNoiseMatrix(state.getDate().durationFrom(prevDate))).getData();
	    estOutput.add(odout);

	    if (odObs.rawMeas.length > 0)
		prevDate = state.getDate();
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
			    if (bias.equalsIgnoreCase(odCfg.parameters.get(k).name))
				biasPos.put(bias, k + 6);
			}
		    }
		}
	    }

	    int Rsize = 0;
	    for (String s: measNames)
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
	    RealMatrix P = new DiagonalMatrix(odCfg.estmCovariance);
	    final Array2DRowRealMatrix sigma = new Array2DRowRealMatrix(numStates, numSigmas);
	    final Array2DRowRealMatrix propSigma = new Array2DRowRealMatrix(numStates, numSigmas);
	    final Array2DRowRealMatrix estimMeas = new Array2DRowRealMatrix(Rsize, numSigmas);
	    final RealMatrix psdCorr = MatrixUtils.createRealIdentityMatrix(P.getRowDimension()).scalarMultiply(1.0E-6);
	    RealMatrix Pprop = null;
	    ArrayRealVector xhat = new ArrayRealVector(xInitial);
	    RealVector xhatPrev = new ArrayRealVector(xInitial);

	    final double[] noiseMult = new double[Rsize];
	    Arrays.fill(noiseMult, 1.0);

	    for (int measIndex = 0; measIndex <= odObs.rawMeas.length; measIndex++)
	    {
		final AbsoluteDate t0 = tm;
		if (measIndex < odObs.rawMeas.length)
		    tm = DataManager.parseDateTime(odObs.rawMeas[measIndex].time);
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

		while (true)
		{
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
			propagator.propagate(stepStart, sigma, stepStart + step, propSigma, enableDMC);
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

		    if (measIndex == odObs.rawMeas.length ||
			(odCfg.propStep != 0.0 && stepStart + step >= bound0 && stepStart + step <= bound1))
		    {
			EstimationOutput odout = new EstimationOutput();
			odout.time = DataManager.getUTCString(new AbsoluteDate(epoch, stepStart + step));
			odout.estimatedState = xhatPrev.toArray();
			odout.propagatedCovariance = Pprop.getData();
			estOutput.add(odout);
		    }

		    stepStart += step;
		    if (FastMath.abs(step) < 1.0E-6 || FastMath.abs(stepFinal - stepStart) < 1.0E-6)
			break;
		}

		if (measIndex == odObs.rawMeas.length)
		    break;

		enableDMC = true;
		RealVector rawMeas = null;
		for (int i = 0; i < numSigmas; i++)
		{
		    double[] pv = propSigma.getColumn(i);
		    ssta[0] = new SpacecraftState(new CartesianOrbit(new PVCoordinates(new Vector3D(pv[0], pv[1], pv[2]),
										       new Vector3D(pv[3], pv[4], pv[5])),
								     odCfg.propFrame, tm, Constants.EGM96_EARTH_MU),
						  propagator.getAttitude(tm, pv), odCfg.rsoMass);

		    if (combinedMeas || measNames.length == 1)
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

		if (odObs.rawMeas[measIndex].station != null)
		{
		    String name = new StringBuilder(odObs.rawMeas[measIndex].station).append(measNames[0]).toString();
		    Integer pos = biasPos.get(name);
		    if (pos != null)
		    {
			if (measNames.length == 2)
			{
			    name = new StringBuilder(odObs.rawMeas[measIndex].station).append(measNames[1]).toString();
			    rawMeas = rawMeas.subtract(new ArrayRealVector(new double[]{xhatPrev.getEntry(pos),
											xhatPrev.getEntry(biasPos.get(name))}));
			}
			else
			    rawMeas = rawMeas.subtract(new ArrayRealVector(new double[]{xhatPrev.getEntry(pos)}));
		    }
		}

		for (int attempt = 0; attempt < 2; attempt++)
		{
		    RealMatrix Pyy = R.copy();
		    for (int i = 0; attempt == 1 && i < Rsize; i++)
			Pyy.setEntry(i, i, Pyy.getEntry(i, i)*noiseMult[i]);

		    RealMatrix Pxy = new Array2DRowRealMatrix(numStates, Rsize);
		    RealVector yhatpre = addColumns(estimMeas).mapMultiplyToSelf(weight);
		    for (int i = 0; i < numSigmas; i++)
		    {
			RealVector y = estimMeas.getColumnVector(i).subtract(yhatpre);
			Pyy = Pyy.add(y.outerProduct(y).scalarMultiply(weight));
			Pxy = Pxy.add(propSigma.getColumnVector(i).subtract(xhatPrev).outerProduct(y).scalarMultiply(weight));
		    }

		    RealMatrix K = Pxy.multiply(MatrixUtils.inverse(Pyy));
		    xhat = new ArrayRealVector(xhatPrev.add(odCfg.parameterMatrix.multiply(K).operate(rawMeas.subtract(yhatpre))));
		    P = Pprop.subtract(odCfg.parameterMatrix.multiply(K.multiply(Pyy.multiply(K.transpose()))));

		    double[] pv = xhat.toArray();
		    ssta[0] = new SpacecraftState(new CartesianOrbit(new PVCoordinates(new Vector3D(pv[0], pv[1], pv[2]),
										       new Vector3D(pv[3], pv[4], pv[5])),
								     odCfg.propFrame, tm, Constants.EGM96_EARTH_MU),
						  propagator.getAttitude(tm, pv), odCfg.rsoMass);

		    EstimationOutput odout = new EstimationOutput();
		    odout.preFit = new HashMap<String, double[]>();
		    odout.postFit = new HashMap<String, double[]>();

		    if (combinedMeas || measNames.length == 1)
		    {
			for (int i = 0; i < measNames.length; i++)
			{
			    double[] fitv = odObs.measObjs.get(measIndex).estimate(0, 0, ssta).getEstimatedValue();
			    if (measNames.length == 1)
			    {
				odout.preFit.put(measNames[i], yhatpre.toArray());
				odout.postFit.put(measNames[i], fitv);
			    }
			    else
			    {
				odout.preFit.put(measNames[i], new double[] {yhatpre.getEntry(i)});
				odout.postFit.put(measNames[i], new double[] {fitv[i]});
			    }
			}
		    }
		    else
		    {
			double[] fitv = odObs.measObjs.get(measIndex*2).estimate(0, 0, ssta).getEstimatedValue();
			odout.preFit.put(measNames[0], new double[] {yhatpre.getEntry(0)});
			odout.postFit.put(measNames[0], fitv);
			fitv = odObs.measObjs.get(measIndex*2 + 1).estimate(0, 0, ssta).getEstimatedValue();
			odout.preFit.put(measNames[1], new double[] {yhatpre.getEntry(1)});
			odout.postFit.put(measNames[1], fitv);
		    }

		    if (attempt == 0 && odCfg.estmOutlierWarmup > 0 && measIndex >= odCfg.estmOutlierWarmup &&
			odCfg.estmOutlierSigma > 0.0)
		    {
			int pos = 0;
			boolean isOutlier = false;
			for (int i = 0; i < measNames.length; i++)
			{
			    double[] fit = odout.postFit.get(measNames[i]);
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

		    odout.time = odObs.rawMeas[measIndex].time;
		    odout.station = odObs.rawMeas[measIndex].station;
		    odout.estimatedState = pv;
		    odout.propagatedCovariance = Pprop.getData();
		    odout.innovationCovariance = Pyy.getData();
		    odout.estimatedCovariance = P.getData();
		    estOutput.add(odout);
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
