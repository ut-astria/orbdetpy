/*
 * Estimation.java - Implementation of estimation algorithms.
 * Copyright (C) 2018 Shiva Iyer <shiva.iyer AT utexas DOT edu>
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

import com.google.gson.GsonBuilder;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import org.astria.DataManager;
import org.astria.ManualPropagation;
import org.astria.Measurements;
import org.astria.PropagatorBuilder;
import org.astria.Settings;
import org.hipparchus.geometry.euclidean.threed.Vector3D;
import org.hipparchus.linear.Array2DRowRealMatrix;
import org.hipparchus.linear.ArrayRealVector;
import org.hipparchus.linear.CholeskyDecomposition;
import org.hipparchus.linear.DiagonalMatrix;
import org.hipparchus.linear.MatrixUtils;
import org.hipparchus.linear.RealMatrix;
import org.hipparchus.linear.RealVector;
import org.orekit.estimation.measurements.GroundStation;
import org.orekit.estimation.measurements.ObservedMeasurement;
import org.orekit.estimation.measurements.Range;
import org.orekit.estimation.sequential.ConstantProcessNoise;
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
import org.orekit.propagation.numerical.NumericalPropagator;
import org.orekit.time.AbsoluteDate;
import org.orekit.time.DateTimeComponents;
import org.orekit.utils.Constants;
import org.orekit.utils.ParameterDriver;
import org.orekit.utils.ParameterDriversList;
import org.orekit.utils.PVCoordinates;
import org.orekit.utils.TimeStampedPVCoordinates;

public class Estimation
{
    protected Settings odcfg;
    protected Measurements odobs;

    protected String[] meanames;
    protected boolean angles;

    protected JSONResults results;

    public final static String DMC_ACC_ESTM = "DMCaccest";
    public final static String DMC_ACC_PROP = "DMCaccprop";

    public Estimation(String cfgjson, String obsjson) throws Exception
    {
	odcfg = Settings.loadJSON(cfgjson);

	if (odcfg.Estimation.Filter == null)
	    throw(new Exception("Missing configuration parameter Estimation.Filter"));
	if (odcfg.Estimation.NoiseTimeDelta <= 0.0)
	    throw(new Exception("Invalid configuration parameter Estimation.NoiseTimeDelta"));

	if (odcfg.Estimation.Filter.equals("UKF"))
	    odcfg.forces.add(0, new NewtonianAttraction(Constants.EGM96_EARTH_MU));

	odobs = Measurements.loadJSON(obsjson, odcfg.stations, odcfg.Measurements);

	meanames = odcfg.Measurements.keySet().toArray(new String[0]);
	angles = meanames[0].equals("Azimuth") || meanames[0].equals("Elevation") ||
	    meanames[0].equals("RightAscension") || meanames[0].equals("Declination");
    }

    public String determineOrbit() throws Exception
    {
	results = new JSONResults();
	results.Filter = odcfg.Estimation.Filter;

	if (odcfg.Estimation.Filter.equals("UKF"))
	    new UnscentedKalmanFilter().determineOrbit();
	else
	    new ExtendedKalmanFilter().determineOrbit();

	return(new GsonBuilder().setPrettyPrinting().create().toJson(results));
    }

    protected class ExtendedKalmanFilter implements KalmanObserver
    {
	protected void determineOrbit()
	{
	    double[] Xi = odcfg.getInitialState();
	    CartesianOrbit X0 = new CartesianOrbit(new PVCoordinates(
						       new Vector3D(Xi[0], Xi[1], Xi[2]),
						       new Vector3D(Xi[3], Xi[4], Xi[5])),
						   DataManager.eme2000, new AbsoluteDate(
						       DateTimeComponents.parseDateTime(odcfg.Propagation.Start),
						       DataManager.utcscale), Constants.EGM96_EARTH_MU);

	    PropagatorBuilder prop = new PropagatorBuilder(odcfg, X0, new DormandPrince853IntegratorBuilder(
							       1E-3, 300.0, 1.0), PositionAngle.MEAN, 10.0);
	    prop.setMass(odcfg.SpaceObject.Mass);
	    for (ForceModel fm : odcfg.forces)
		prop.addForceModel(fm);

	    ParameterDriversList plst = prop.getPropagationParametersDrivers();
	    for (Settings.EstimatedParameter ep : odcfg.estparams)
	    {
		ParameterDriver pdrv = new ParameterDriver(ep.name, ep.value, 1.0, ep.min, ep.max);
		pdrv.setSelected(true);
		plst.add(pdrv);
	    }

	    KalmanEstimatorBuilder build = new KalmanEstimatorBuilder();
	    build.addPropagationConfiguration(prop, new ConstantProcessNoise(new DiagonalMatrix(odcfg.Estimation.Covariance),
									     odcfg.getProcessNoiseMatrix()));

	    KalmanEstimator filt = build.build();
	    filt.setObserver(this);
	    NumericalPropagator est = filt.processMeasurements(odobs.measobjs)[0];

	    results.Propagation.Time = odcfg.Propagation.End;
	    results.Propagation.State = new double[odcfg.estparams.size() + 6];
	    PVCoordinates pvc = est.getPVCoordinates(new AbsoluteDate(DateTimeComponents.parseDateTime(odcfg.Propagation.End),
								      DataManager.utcscale), DataManager.eme2000);
	    System.arraycopy(pvc.getPosition().toArray(), 0, results.Propagation.State, 0, 3);
	    System.arraycopy(pvc.getVelocity().toArray(), 0, results.Propagation.State, 3, 3);
	    if (odcfg.estparams.size() > 0)
		System.arraycopy(results.Estimation.get(results.Estimation.size() - 1).EstimatedState, 6,
		results.Propagation.State, 6, odcfg.estparams.size());
	}

	public void evaluationPerformed​(KalmanEstimation est)
	{
	    int n = est.getCurrentMeasurementNumber()-1;
	    if (!angles)
		n /= meanames.length;

	    String k;
	    JSONResults.JSONEstimation res;
	    if (results.Estimation.size() <= n)
	    {
		k = meanames[0];
		res = results.new JSONEstimation();
		res.Time = odobs.rawmeas[n].Time;
		results.Estimation.add(res);
	    }
	    else
	    {
		k = meanames[1];
		res = results.Estimation.get(n);
	    }

	    SpacecraftState ssta = est.getPredictedSpacecraftStates()[0];
	    PVCoordinates pvc = ssta.getPVCoordinates();

	    res.EstimatedState = new double[odcfg.estparams.size() + 6];
	    System.arraycopy(pvc.getPosition().toArray(), 0, res.EstimatedState, 0, 3);
	    System.arraycopy(pvc.getVelocity().toArray(), 0, res.EstimatedState, 3, 3);

	    res.EstimatedAcceleration = new double[3];
	    System.arraycopy(pvc.getAcceleration().toArray(), 0, res.EstimatedAcceleration, 0, 3);

	    int i = 6;
	    List<ParameterDriversList.DelegatingDriver> plst = est.getEstimatedPropagationParameters().getDrivers();
	    for (Settings.EstimatedParameter ep : odcfg.estparams)
		for (ParameterDriversList.DelegatingDriver dd : plst)
		    if (dd.getName().equals(ep.name))
			res.EstimatedState[i++] = dd.getValue();

	    if (ssta.hasAdditionalState(Estimation.DMC_ACC_PROP))
	    {
		double[] accric = ssta.getAdditionalState(Estimation.DMC_ACC_PROP);
		System.arraycopy(accric, 0, res.EstimatedState, odcfg.estparams.size() + 3, 3);
	    }

	    double[] pre = est.getPredictedMeasurement().getEstimatedValue();
	    double[] pos = est.getCorrectedMeasurement().getEstimatedValue();
	    if (angles)
	    {
		res.PreFit.put(meanames[0], pre[0]);
		res.PreFit.put(meanames[1], pre[1]);
		res.PostFit.put(meanames[0], pos[0]);
		res.PostFit.put(meanames[1], pos[1]);
		res.InnovationCovariance = est.getPhysicalInnovationCovarianceMatrix().getData();
	    }
	    else
	    {
		res.PreFit.put(k, pre[0]);
		res.PostFit.put(k, pos[0]);
		if (res.InnovationCovariance == null)
		    res.InnovationCovariance = new double[][]{{
			    est.getPhysicalInnovationCovarianceMatrix().getData()[0][0], 0.0}, {0.0, 0.0}};
		else
		    res.InnovationCovariance[1][1] = est.getPhysicalInnovationCovarianceMatrix().getData()[0][0];
	    }
	    res.EstimatedCovariance = est.getPhysicalEstimatedCovarianceMatrix().getData();
	}
    }

    protected class UnscentedKalmanFilter
    {
	protected void determineOrbit() throws Exception
	{
	    int numsta = odcfg.estparams.size() + 6;
	    int numsig = 2*numsta;
	    int veclen = numsta*numsig;

	    ManualPropagation prop = new ManualPropagation(odcfg, veclen);

	    RealMatrix P = new DiagonalMatrix(odcfg.Estimation.Covariance);
	    RealMatrix Q = odcfg.getProcessNoiseMatrix();
	    Array2DRowRealMatrix R = new Array2DRowRealMatrix(odcfg.Measurements.size(), odcfg.Measurements.size());
	    int i = 0, j = 0;
	    for (i = 0; i < 2; i++)
	    {
		Settings.JSONMeasurement m = odcfg.Measurements.get(meanames[i]);
		R.setEntry(i, i, m.Error*m.Error);
	    }

	    double[] Xi = odcfg.getInitialState();
	    Array2DRowRealMatrix sigma = new Array2DRowRealMatrix(numsta, numsig);
	    Array2DRowRealMatrix sigpr = new Array2DRowRealMatrix(numsta, numsig);
	    Array2DRowRealMatrix spupd = new Array2DRowRealMatrix(odcfg.Measurements.size(), numsig);
	    ArrayRealVector xhat = new ArrayRealVector(Xi);
	    RealVector xhatpre = new ArrayRealVector(Xi);
	    double weight = 0.5/numsta;
	    double[] spvec = new double[veclen];

	    double tofm = 0.0;
	    AbsoluteDate epoch = new AbsoluteDate(DateTimeComponents.parseDateTime(odcfg.Propagation.Start),
						  DataManager.utcscale);
	    AbsoluteDate tm = new AbsoluteDate(epoch, 0.0);
	    SpacecraftState[] ssta = new SpacecraftState[1];

	    for (int mix = 0; mix <= odobs.rawmeas.length; mix++)
	    {
		JSONResults.JSONEstimation odout = results.new JSONEstimation();

		double tof0 = tofm;
		AbsoluteDate t0 = new AbsoluteDate(tm, 0.0);
		if (mix < odobs.rawmeas.length)
		{
		    tm = new AbsoluteDate(DateTimeComponents.parseDateTime(odobs.rawmeas[mix].Time),
					  DataManager.utcscale);
		    AbsoluteDate tmlt = new AbsoluteDate(tm, -tofm);

		    double[] pv = xhat.toArray();
		    TimeStampedPVCoordinates pvs = new TimeStampedPVCoordinates(tmlt,
										new Vector3D(pv[0], pv[1], pv[2]),
										new Vector3D(pv[3], pv[4], pv[5]));
		    PVCoordinates pvi = odcfg.stations.get(odobs.rawmeas[mix].Station).getBaseFrame().
			getPVCoordinates(tm, DataManager.eme2000);
		    tofm = Range.signalTimeOfFlight(pvs, pvi.getPosition(), tm);

		    odout.Time = odobs.rawmeas[mix].Time;
		}
		else
		{
		    tofm = 0.0;
		    tm = new AbsoluteDate(DateTimeComponents.parseDateTime(odcfg.Propagation.End),
					  DataManager.utcscale);
		}

		RealMatrix sqrP = new CholeskyDecomposition(P.scalarMultiply(numsta), 1E-6,
							    CholeskyDecomposition.DEFAULT_ABSOLUTE_POSITIVITY_THRESHOLD).getL();
		for (i = 0; i < numsta; i++)
		{
		    sigma.setColumnVector(i, xhat.add(sqrP.getColumnVector(i)));
		    sigma.setColumnVector(numsta+i, xhat.subtract(sqrP.getColumnVector(i)));
		}

		for (j = 0; j < odcfg.estparams.size(); j++)
		{
		    Settings.EstimatedParameter tempep = odcfg.estparams.get(j);
		    for (i = 0; i < numsig; i++)
		    {
			double tempdbl = sigma.getEntry(j + 6, i);
			sigma.setEntry(j + 6, i, Math.min(Math.max(tempdbl, tempep.min), tempep.max));
		    }
		}

		unstack(sigpr, prop.propagate(t0.durationFrom(epoch) - tof0, stack(sigma, spvec),
					      tm.durationFrom(epoch) - tofm));
		xhatpre = addColumns(sigpr).mapMultiplyToSelf(weight);
		if (mix == odobs.rawmeas.length)
		    break;

		RealVector raw = null;
		RealMatrix Ppre = Q.copy();
		AbsoluteDate tmlt = new AbsoluteDate(tm, -tofm);
		for (i = 0; i < numsig; i++)
		{
		    RealVector y = sigpr.getColumnVector(i).subtract(xhatpre);
		    Ppre = Ppre.add(y.outerProduct(y).scalarMultiply(weight));

		    double[] pv = sigpr.getColumn(i);
		    ssta[0] = new SpacecraftState(new CartesianOrbit(new PVCoordinates(new Vector3D(pv[0], pv[1], pv[2]),
										       new Vector3D(pv[3], pv[4], pv[5])),
								     DataManager.eme2000, tmlt, Constants.EGM96_EARTH_MU),
						  odcfg.SpaceObject.Mass);

		    if (angles)
		    {
			double[] fitv = odobs.measobjs.get(mix).estimate(1, 1, ssta).getEstimatedValue();
			spupd.setColumn(i, fitv);
			if (raw == null)
			    raw = new ArrayRealVector(odobs.measobjs.get(mix).getObservedValue());
		    }
		    else
		    {
			double[] fitv = odobs.measobjs.get(mix*2).estimate(1, 1, ssta).getEstimatedValue();
			spupd.setEntry(0, i, fitv[0]);
			fitv = odobs.measobjs.get(mix*2 + 1).estimate(1, 1, ssta).getEstimatedValue();
			spupd.setEntry(1, i, fitv[0]);
			if (raw == null)
			    raw = new ArrayRealVector(new double[]{odobs.measobjs.get(mix*2).getObservedValue()[0],
								   odobs.measobjs.get(mix*2 + 1).getObservedValue()[0]});
		    }
		}

		RealMatrix Pyy = R.copy();
		RealMatrix Pxy = new Array2DRowRealMatrix(numsta, odcfg.Measurements.size());
		RealVector yhatpre = addColumns(spupd).mapMultiplyToSelf(weight);
		for (i = 0; i < numsig; i++)
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
								 DataManager.eme2000, tmlt, Constants.EGM96_EARTH_MU),
						  odcfg.SpaceObject.Mass);

		odout.EstimatedState = pv;
		odout.EstimatedCovariance = P.getData();
		odout.InnovationCovariance = Pyy.getData();
		odout.PreFit.put(meanames[0], yhatpre.getEntry(0));
		odout.PreFit.put(meanames[1], yhatpre.getEntry(1));
		if (angles)
		{
		    double[] fitv = odobs.measobjs.get(mix).estimate(1, 1, ssta).getEstimatedValue();
		    odout.PostFit.put(meanames[0], fitv[0]);
		    odout.PostFit.put(meanames[1], fitv[1]);
		}
		else
		{
		    double[] fitv = odobs.measobjs.get(mix*2).estimate(1, 1, ssta).getEstimatedValue();
		    odout.PostFit.put(meanames[0], fitv[0]);
		    fitv = odobs.measobjs.get(mix*2 + 1).estimate(1, 1, ssta).getEstimatedValue();
		    odout.PostFit.put(meanames[1], fitv[0]);
		}

		results.Estimation.add(odout);
	    }

	    results.Propagation.Time = odcfg.Propagation.End;
	    results.Propagation.State = xhatpre.toArray();
	}

	double[] stack(RealMatrix mat, double[] arr)
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

	RealMatrix unstack(RealMatrix mat, double[] arr)
	{
	    int i,j;
	    int m = mat.getRowDimension();
	    int n = mat.getColumnDimension();

	    for (i = 0; i < n; i++)
		for (j = 0; j < m; j++)
		    mat.setEntry(j, i, arr[i*m + j]);

	    return(mat);
	}

	ArrayRealVector addColumns(RealMatrix mat)
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
    }

    class JSONResults
    {
	class JSONEstimation
	{
	    String Time;
	    HashMap<String, Double> PreFit;
	    HashMap<String, Double> PostFit;
	    double[] EstimatedState;
	    double[] EstimatedAcceleration;
	    double[][] EstimatedCovariance;
	    double[][] InnovationCovariance;

	    public JSONEstimation()
	    {
		PreFit = new HashMap<String, Double>();
		PostFit = new HashMap<String, Double>();
	    }
	}

	class JSONPropagation
	{
	    String Time;
	    double[] State;
	}

	String Filter;
	ArrayList<JSONEstimation> Estimation;
	JSONPropagation Propagation;

	public JSONResults()
	{
	    Estimation = new ArrayList<JSONEstimation>();
	    Propagation = new JSONPropagation();
	}
    }
}
