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
import org.hipparchus.ode.sampling.ODEFixedStepHandler;
import org.hipparchus.ode.sampling.StepNormalizer;
import org.hipparchus.ode.sampling.StepNormalizerBounds;
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
import org.orekit.time.DateTimeComponents;
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

	if (odCfg.estmFilter.equals("UKF") && odCfg.gravityDegree >= 2 && odCfg.gravityOrder >= 0)
	    odCfg.forces.add(0, new NewtonianAttraction(Constants.EGM96_EARTH_MU));

	measNames = odCfg.cfgMeasurements.keySet().toArray(new String[0]);
	Arrays.sort(measNames);
	if (measNames[0].equals("Declination"))
	    measNames = new String[]{"RightAscension", "Declination"};
	combinedMeas = !measNames[0].startsWith("Range");

	epoch = new AbsoluteDate(DateTimeComponents.parseDateTime(odCfg.propStart), DataManager.getTimeScale("UTC"));
	propEnd = new AbsoluteDate(DateTimeComponents.parseDateTime(odCfg.propEnd), DataManager.getTimeScale("UTC"));

	if (odCfg.propStepHandlerStartTime != null)
	    stepHandlerStart = new AbsoluteDate(DateTimeComponents.parseDateTime(odCfg.propStepHandlerStartTime),
						DataManager.getTimeScale("UTC"));
	else
	    stepHandlerStart = epoch;

	if (odCfg.propStepHandlerEndTime != null)
	    stepHandlerEnd = new AbsoluteDate(DateTimeComponents.parseDateTime(odCfg.propStepHandlerEndTime),
					      DataManager.getTimeScale("UTC"));
	else
	    stepHandlerEnd = propEnd;
    }

    public ArrayList<EstimationOutput> determineOrbit()
    {
	int size = FastMath.max(odObs.rawMeas.length, 10);
	if (odCfg.propStep != 0.0)
	    size += (int) FastMath.abs(propEnd.durationFrom(epoch)/odCfg.propStep) + 2;
	estOutput = new ArrayList<EstimationOutput>(size);

	if (odCfg.estmFilter.equals("UKF"))
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
	    double[] Xi = odCfg.getInitialState();
	    prevDate = epoch;
	    prevPosition = new Vector3D(Xi[0],Xi[1],Xi[2]);
	    prevVelocity = new Vector3D(Xi[3],Xi[4],Xi[5]);
	    prevCovariance = new DiagonalMatrix(odCfg.estmCovariance);
	    CartesianOrbit X0 = new CartesianOrbit(new PVCoordinates(prevPosition, prevVelocity),
						   odCfg.propFrame, epoch, Constants.EGM96_EARTH_MU);

	    OrekitFixedStepHandler handler = null;
    	    if (odCfg.propStep > 0.0)
		handler = this;

	    propBuilder = new PropagatorBuilder(odCfg, X0, new DormandPrince853IntegratorBuilder(
						    odCfg.integMinTimeStep, odCfg.integMaxTimeStep, 1.0), PositionAngle.TRUE, 10.0, handler);
	    propBuilder.setMass(odCfg.rsoMass);
	    for (ForceModel fm : odCfg.forces)
		propBuilder.addForceModel(fm);

	    ParameterDriversList plst = propBuilder.getPropagationParametersDrivers();
	    for (Settings.Parameter ep : odCfg.parameters)
	    {
		ParameterDriver pdrv = new ParameterDriver(ep.name, ep.value, 1.0, ep.min, ep.max);
		pdrv.setSelected(true);
		plst.add(pdrv);
	    }

	    AttitudeProvider attprov = odCfg.getAttitudeProvider();
	    if (attprov != null)
		propBuilder.setAttitudeProvider(attprov);

	    KalmanEstimatorBuilder builder = new KalmanEstimatorBuilder();
	    builder.addPropagationConfiguration(propBuilder, this);
	    KalmanEstimator filter = builder.build();
	    filter.setObserver(this);

	    AbstractIntegratedPropagator estimator = null;
	    AbstractIntegratedPropagator[] estimators = filter.processMeasurements(odObs.measObjs);
	    if (estimators != null)
		estimator = estimators[0];
	    else
		estimator = propBuilder.buildPropagator(propBuilder.getSelectedNormalizedParameters());

	    if (odObs.rawMeas.length == 0 || !odCfg.propEnd.equals(odObs.rawMeas[odObs.rawMeas.length-1].time))
	    {
		SpacecraftState state = estimator.propagate(propEnd);
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
	    double tmeas = FastMath.abs(curr.getDate().durationFrom(prev.getDate()));
	    if (tmeas > 0.0)
		measDeltaT = tmeas;
	    else
		tmeas = measDeltaT;
	    return(odCfg.getProcessNoiseMatrix(tmeas));
	}

	@Override public void evaluationPerformed(KalmanEstimation est)
	{
	    int n = est.getCurrentMeasurementNumber() - 1;
	    if (!combinedMeas)
		n /= measNames.length;
	    Measurements.Measurement raw = odObs.rawMeas[n];

	    EstimationOutput res = null;
	    for (EstimationOutput loop: estOutput)
	    {
		if (loop.preFit != null && loop.postFit != null && loop.time.equals(raw.time) &&
		    (loop.station == null || loop.station.equals(raw.station)))
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

	    SpacecraftState ssta = est.getPredictedSpacecraftStates()[0];
	    PVCoordinates pvc = ssta.getPVCoordinates();
	    res.estimatedState = new double[odCfg.parameters.size() + 6];
	    System.arraycopy(pvc.getPosition().toArray(), 0, res.estimatedState, 0, 3);
	    System.arraycopy(pvc.getVelocity().toArray(), 0, res.estimatedState, 3, 3);

	    ParameterDriversList plst = est.getEstimatedPropagationParameters();
	    for (int i = 0; i < odCfg.parameters.size(); i++)
	    {
		Settings.Parameter ep = odCfg.parameters.get(i);
		res.estimatedState[i + 6] = est.getEstimatedPropagationParameters().findByName(ep.name).getValue();
	    }

	    double[] pre = est.getPredictedMeasurement().getEstimatedValue();
	    double[] pos = est.getCorrectedMeasurement().getEstimatedValue();
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
		res.innovationCovariance = est.getPhysicalInnovationCovarianceMatrix().getData();
	    }
	    else
	    {
		res.preFit.put(key, pre);
		res.postFit.put(key, pos);
		if (res.innovationCovariance == null)
		    res.innovationCovariance = new double[][]{{
			    est.getPhysicalInnovationCovarianceMatrix().getData()[0][0], 0.0}, {0.0, 0.0}};
		else
		    res.innovationCovariance[1][1] = est.getPhysicalInnovationCovarianceMatrix().getData()[0][0];
	    }

	    RealMatrix phi = est.getPhysicalStateTransitionMatrix();
	    res.propagatedCovariance = phi.multiply(prevCovariance).multiply(phi.transpose()).add(
		odCfg.getProcessNoiseMatrix(measDeltaT)).getData();

	    res.estimatedCovariance = est.getPhysicalEstimatedCovarianceMatrix().getData();
	    prevCovariance = est.getPhysicalEstimatedCovarianceMatrix();
	}

	@Override public void handleStep(SpacecraftState state, boolean lastStep)
	{
	    if (state.getDate().durationFrom(stepHandlerStart) < 0.0 || state.getDate().durationFrom(stepHandlerEnd) > 0.0)
		return;

	    PVCoordinates pvc = state.getPVCoordinates(odCfg.propFrame);
	    for (ObservedMeasurement m: odObs.measObjs)
	    {
		if (m.getDate().durationFrom(state.getDate()) == 0.0)
		    return;
	    }

	    EstimationOutput odout = new EstimationOutput();
	    odout.time = DataManager.getUTCString(state.getDate());
	    odout.estimatedState = new double[odCfg.parameters.size() + 6];
	    System.arraycopy(pvc.getPosition().toArray(), 0, odout.estimatedState, 0, 3);
	    System.arraycopy(pvc.getVelocity().toArray(), 0, odout.estimatedState, 3, 3);
	    if (odCfg.parameters.size() > 0 && estOutput.size() > 0)
		System.arraycopy(estOutput.get(estOutput.size() - 1).estimatedState, 6,
				 odout.estimatedState, 6, odCfg.parameters.size());

	    Rotation phi1 = new Rotation(prevPosition, pvc.getPosition());
	    Rotation phi2 = new Rotation(prevVelocity, pvc.getVelocity());
	    RealMatrix phi = MatrixUtils.createRealIdentityMatrix(prevCovariance.getRowDimension());
	    phi.setSubMatrix(phi1.getMatrix(), 0, 0);
	    phi.setSubMatrix(phi2.getMatrix(), 3, 3);
	    prevPosition = pvc.getPosition();
	    prevVelocity = pvc.getVelocity();

	    odout.propagatedCovariance = phi.multiply(prevCovariance).multiply(phi.transpose()).add(
		odCfg.getProcessNoiseMatrix(FastMath.abs(state.getDate().durationFrom(prevDate)))).getData();
	    estOutput.add(odout);
	    prevDate = state.getDate();
	}
    }

    private class UnscentedKalmanFilter
    {
	private void determineOrbit()
	{
	    int numsta = odCfg.parameters.size() + 6;
	    int numsig = 2*numsta;
	    int veclen = numsta*numsig;
	    double weight = 0.5/numsta;
	    RealMatrix P = new DiagonalMatrix(odCfg.estmCovariance);
	    double[] Xi = odCfg.getInitialState();

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

	    StepNormalizer normalizer = null;
	    IntegrationStepHandler handler = new IntegrationStepHandler();
    	    if (odCfg.propStep > 0.0)
		normalizer = new StepNormalizer(odCfg.propStep, handler, StepNormalizerBounds.BOTH);

	    Array2DRowRealMatrix sigma = new Array2DRowRealMatrix(numsta, numsig);
	    Array2DRowRealMatrix sigpr = new Array2DRowRealMatrix(numsta, numsig);
	    Array2DRowRealMatrix spupd = new Array2DRowRealMatrix(Rsize, numsig);
	    RealMatrix psdCorr = MatrixUtils.createRealIdentityMatrix(P.getRowDimension()).scalarMultiply(1.0E-6);
	    ArrayRealVector xhat = new ArrayRealVector(Xi);
	    RealVector xhatpre = new ArrayRealVector(Xi);
	    double[] spvec = new double[veclen];
	    AbsoluteDate tm = new AbsoluteDate(epoch, 0.0);
	    SpacecraftState[] ssta = new SpacecraftState[1];

	    ManualPropagation prop = new ManualPropagation(odCfg, veclen, normalizer);
	    prop.setDMCState(false);

	    for (int mix = 0; mix <= odObs.rawMeas.length; mix++)
	    {
		AbsoluteDate t0 = new AbsoluteDate(tm, 0.0);
		if (mix < odObs.rawMeas.length)
		    tm = new AbsoluteDate(DateTimeComponents.parseDateTime(odObs.rawMeas[mix].time), DataManager.getTimeScale("UTC"));
		else
		    tm = propEnd;

		RealMatrix Ptemp = P.scalarMultiply(numsta);
		RealMatrix sqrP = new CholeskyDecomposition(
		    Ptemp.add(Ptemp.transpose()).scalarMultiply(0.5).add(psdCorr), 1E-6, 1E-16).getL();
		for (int i = 0; i < numsta; i++)
		{
		    sigma.setColumnVector(i, xhat.add(sqrP.getColumnVector(i)));
		    sigma.setColumnVector(numsta + i, xhat.subtract(sqrP.getColumnVector(i)));
		}

		if (odCfg.parameters.size() > 0)
		{
		    double[][] sigdata = sigma.getData();
		    for (int j = 6; j < odCfg.parameters.size() + 6; j++)
		    {
			Settings.Parameter tempep = odCfg.parameters.get(j - 6);
			for (int i = 0; i < numsig; i++)
			    sigdata[j][i] = FastMath.min(FastMath.max(sigdata[j][i], tempep.min), tempep.max);
		    }
		    sigma.setSubMatrix(sigdata, 0, 0);
		}

		double propt0 = t0.durationFrom(epoch);
		double propt1 = tm.durationFrom(epoch);
		if (propt0 == propt1)
		    sigpr.setSubMatrix(sigma.getData(), 0, 0);
		else
		{
		    ODEStateAndDerivative state = prop.propagate(propt0, stackColumns(sigma, spvec), propt1);
		    unstackColumn(sigpr, state.getPrimaryState());
		    prop.setDMCState(true);
		    if (normalizer == null)
			handler.handleStep(state, true);
		}

		if (mix == odObs.rawMeas.length)
		    break;

		RealVector raw = null;
		xhatpre = addColumns(sigpr).mapMultiplyToSelf(weight);
		RealMatrix Ppre = odCfg.getProcessNoiseMatrix(FastMath.abs(propt1 - propt0));
		for (int i = 0; i < numsig; i++)
		{
		    RealVector y = sigpr.getColumnVector(i).subtract(xhatpre);
		    Ppre = Ppre.add(y.outerProduct(y).scalarMultiply(weight));
		    double[] pv = sigpr.getColumn(i);
		    ssta[0] = new SpacecraftState(new CartesianOrbit(new PVCoordinates(new Vector3D(pv[0], pv[1], pv[2]),
										       new Vector3D(pv[3], pv[4], pv[5])),
								     odCfg.propFrame, tm, Constants.EGM96_EARTH_MU),
						  prop.getAttitude(tm, pv), odCfg.rsoMass);

		    if (combinedMeas || measNames.length == 1)
		    {
			double[] fitv = odObs.measObjs.get(mix).estimate(0, 0, ssta).getEstimatedValue();
			spupd.setColumn(i, fitv);
			if (raw == null)
			    raw = new ArrayRealVector(odObs.measObjs.get(mix).getObservedValue());
		    }
		    else
		    {
			double[] fitv = odObs.measObjs.get(mix*2).estimate(0, 0, ssta).getEstimatedValue();
			spupd.setEntry(0, i, fitv[0]);
			fitv = odObs.measObjs.get(mix*2 + 1).estimate(0, 0, ssta).getEstimatedValue();
			spupd.setEntry(1, i, fitv[0]);
			if (raw == null)
			    raw = new ArrayRealVector(new double[]{odObs.measObjs.get(mix*2).getObservedValue()[0],
								   odObs.measObjs.get(mix*2 + 1).getObservedValue()[0]});
		    }
		}

		if (odObs.rawMeas[mix].station != null)
		{
		    String name = new StringBuilder(odObs.rawMeas[mix].station).append(measNames[0]).toString();
		    Integer pos = biasPos.get(name);
		    if (pos != null)
		    {
			if (measNames.length == 2)
			{
			    name = new StringBuilder(odObs.rawMeas[mix].station).append(measNames[1]).toString();
			    raw = raw.subtract(new ArrayRealVector(new double[]{xhatpre.getEntry(pos),
										xhatpre.getEntry(biasPos.get(name))}));
			}
			else
			    raw = raw.subtract(new ArrayRealVector(new double[]{xhatpre.getEntry(pos)}));
		    }
		}

		RealMatrix Pyy = R.copy();
		RealMatrix Pxy = new Array2DRowRealMatrix(numsta, Rsize);
		RealVector yhatpre = addColumns(spupd).mapMultiplyToSelf(weight);
		for (int i = 0; i < numsig; i++)
		{
		    RealVector y = spupd.getColumnVector(i).subtract(yhatpre);
		    Pyy = Pyy.add(y.outerProduct(y).scalarMultiply(weight));
		    Pxy = Pxy.add(sigpr.getColumnVector(i).subtract(xhatpre).outerProduct(y).scalarMultiply(weight));
		}

		RealMatrix K = Pxy.multiply(MatrixUtils.inverse(Pyy));
		xhat = new ArrayRealVector(xhatpre.add(K.operate(raw.subtract(yhatpre))));
		P = Ppre.subtract(K.multiply(Pyy.multiply(K.transpose())));

		double[] pv = xhat.toArray();
		ssta[0] = new SpacecraftState(new CartesianOrbit(new PVCoordinates(new Vector3D(pv[0], pv[1], pv[2]),
										   new Vector3D(pv[3], pv[4], pv[5])),
								 odCfg.propFrame, tm, Constants.EGM96_EARTH_MU),
					      prop.getAttitude(tm, pv), odCfg.rsoMass);

		EstimationOutput odout = new EstimationOutput();
		odout.time = odObs.rawMeas[mix].time;
		odout.station = odObs.rawMeas[mix].station;
		odout.estimatedState = pv;
		odout.propagatedCovariance = Ppre.getData();
		odout.innovationCovariance = Pyy.getData();
		odout.estimatedCovariance = P.getData();
		odout.preFit = new HashMap<String, double[]>();
		odout.postFit = new HashMap<String, double[]>();

		if (combinedMeas || measNames.length == 1)
		{
		    for (int i = 0; i < measNames.length; i++)
		    {
			double[] fitv = odObs.measObjs.get(mix).estimate(0, 0, ssta).getEstimatedValue();
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
		    double[] fitv = odObs.measObjs.get(mix*2).estimate(0, 0, ssta).getEstimatedValue();
		    odout.preFit.put(measNames[0], new double[] {yhatpre.getEntry(0)});
		    odout.postFit.put(measNames[0], fitv);
		    fitv = odObs.measObjs.get(mix*2 + 1).estimate(0, 0, ssta).getEstimatedValue();
		    odout.preFit.put(measNames[1], new double[] {yhatpre.getEntry(1)});
		    odout.postFit.put(measNames[1], fitv);
		}

		estOutput.add(odout);
	    }
	}

	private double[] stackColumns(RealMatrix mat, double[] arr)
	{
	    int i,j;
	    int m = mat.getRowDimension();
	    int n = mat.getColumnDimension();
	    double[][] matdata = mat.getData();
	    for (i = 0; i < n; i++)
		for (j = 0; j < m; j++)
		    arr[i*m + j] = matdata[j][i];

	    return(arr);
	}

	private RealMatrix unstackColumn(RealMatrix mat, double[] arr)
	{
	    int m = mat.getRowDimension();
	    int n = mat.getColumnDimension();
	    for (int i = 0; i < n; i++)
		mat.setColumn(i, Arrays.copyOfRange(arr, i*m, (i+1)*m));

	    return(mat);
	}

	private ArrayRealVector addColumns(RealMatrix mat)
	{
	    int i,j;
	    double sum;
	    int m = mat.getRowDimension();
	    int n = mat.getColumnDimension();
	    double[][] arr = mat.getData();
	    ArrayRealVector out = new ArrayRealVector(m);
	    for (j = 0; j < m; j++)
	    {
		sum = 0.0;
		for (i = 0; i < n; i++)
		    sum += arr[j][i];
		out.setEntry(j, sum);
	    }

	    return(out);
	}

	private class IntegrationStepHandler implements ODEFixedStepHandler
	{
	    private final int numStates;
	    private final double weight;
	    private final Array2DRowRealMatrix propSigma;
	    private double prevTime;

	    public IntegrationStepHandler()
	    {
		numStates = odCfg.parameters.size() + 6;
		weight = 0.5/numStates;
		propSigma = new Array2DRowRealMatrix(numStates, 2*numStates);
	    }

	    @Override public void handleStep(ODEStateAndDerivative state, boolean lastStep)
	    {
		AbsoluteDate time = new AbsoluteDate(epoch, state.getTime());
		if (time.durationFrom(stepHandlerStart) < 0.0 || time.durationFrom(stepHandlerEnd) > 0.0)
		    return;

		for (ObservedMeasurement m: odObs.measObjs)
		{
		    if (m.getDate().durationFrom(time) == 0.0)
			return;
		}

		EstimationOutput odout = new EstimationOutput();
		odout.time = DataManager.getUTCString(time);
		odout.estimatedState = new double[numStates];

		double[] X = state.getPrimaryState();
		unstackColumn(propSigma, X);
		RealVector xhatpre = addColumns(propSigma).mapMultiplyToSelf(weight);
		System.arraycopy(xhatpre.toArray(), 0, odout.estimatedState, 0, odout.estimatedState.length);

		RealMatrix Ppre = odCfg.getProcessNoiseMatrix(FastMath.abs(state.getTime() - prevTime));
		for (int i = 0; i < 2*numStates; i++)
		{
		    RealVector y = propSigma.getColumnVector(i).subtract(xhatpre);
		    Ppre = Ppre.add(y.outerProduct(y).scalarMultiply(weight));
		}

		odout.propagatedCovariance = Ppre.getData();
		estOutput.add(odout);
		prevTime = state.getTime();
	    }
	}
    }
}
