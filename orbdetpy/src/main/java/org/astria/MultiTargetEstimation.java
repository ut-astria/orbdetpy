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

public final class MultiTargetEstimation
{
    private Settings odCfg;

    private final Measurements odObs;

    private String[] measNames;
    private final boolean combinedMeas;

    private final AbsoluteDate epoch;
    private final AbsoluteDate propEnd;
    private final AbsoluteDate stepHandlerStart;
    private final AbsoluteDate stepHandlerEnd;
    private ArrayList<Estimation.EstimationOutput> estOutput;

    public MultiTargetEstimation(ArrayList<Settings> cfgList, ArrayList<Measurements> obsList)
    {
	this(cfgList.get(0), obsList.get(0));
    }

    public MultiTargetEstimation(Settings odCfg, Measurements odObs)
    {
	odCfg.forces.add(0, new NewtonianAttraction(Constants.EGM96_EARTH_MU));
    	
	this.odCfg = odCfg;
	this.odObs = odObs;
	
	measNames = this.odCfg.cfgMeasurements.keySet().toArray(new String[0]);
	Arrays.sort(measNames);
	if (measNames[0].equalsIgnoreCase("Declination"))
	    measNames = new String[]{"RightAscension", "Declination"};
	combinedMeas = !measNames[0].equalsIgnoreCase("Range") && !measNames[0].equalsIgnoreCase("RangeRate");

	epoch = DataManager.parseDateTime(this.odCfg.propStart);
	propEnd = DataManager.parseDateTime(this.odCfg.propEnd);

	if (this.odCfg.propStepHandlerStartTime != null)
	    stepHandlerStart = DataManager.parseDateTime(this.odCfg.propStepHandlerStartTime);
	else
	{
	    if (this.odCfg.propStep > 0.0)
		stepHandlerStart = epoch;
	    else
		stepHandlerStart = propEnd;
	}

	if (this.odCfg.propStepHandlerEndTime != null)
	    stepHandlerEnd = DataManager.parseDateTime(this.odCfg.propStepHandlerEndTime);
	else
	{
	    if (this.odCfg.propStep > 0.0)
		stepHandlerEnd = propEnd;
	    else
		stepHandlerEnd = epoch;
	}
    }

    public ArrayList<Estimation.EstimationOutput> multiTargetDetermineOrbit()
    {
	int size = FastMath.max(odObs.rawMeas.length, 10);
	if (this.odCfg.propStep != 0.0)
	    size += (int) FastMath.abs(propEnd.durationFrom(epoch)/this.odCfg.propStep) + 2;
	estOutput = new ArrayList<Estimation.EstimationOutput>(size);

    new MultiTargetFilter().determineOrbit();

    
	return(estOutput);
    }

    private class MultiTargetFilter
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
	    
	    final int totalObjNum = 1;
	    final int numStates = odCfg.parameters.size() + 6;
	    final int numSigmas = 2*numStates;
	    final double weight = 0.5/numStates;
	    final double bound0 = stepHandlerStart.durationFrom(epoch);
	    final double bound1 = stepHandlerEnd.durationFrom(epoch);
	    AbsoluteDate tm = epoch;
	    final SpacecraftState[] ssta = new SpacecraftState[1];
	    final ManualPropagation propagator = new ManualPropagation(odCfg);

	    
	    ArrayList<SingleObject> Objects = new ArrayList<SingleObject>();
	    
	    //Create object for each RSO
		for(int objNum = 0; objNum < totalObjNum; objNum++)
		{
			Objects.add(new SingleObject(numStates, numSigmas, Rsize));
			
		}
		
	    
		EstimationLoop:
	    for (int measIndex = 0; measIndex <= odObs.rawMeas.length; measIndex++)
	    {
		final AbsoluteDate t0 = tm;
		if (measIndex < odObs.rawMeas.length)
		    tm = DataManager.parseDateTime(odObs.rawMeas[measIndex].time);
		else
		    tm = propEnd;

		double propStart = t0.durationFrom(epoch);
		final double propFinal = tm.durationFrom(epoch);
		double stepSum = 0.0;
		
		//propagate
		for(int objNum = 0; objNum < totalObjNum; objNum++)
		{
			
			SingleObject currSC = Objects.get(objNum);
			
			final RealMatrix Ptemp = currSC.P.scalarMultiply(numStates);
			final RealMatrix sqrtP = new CholeskyDecomposition(
			    Ptemp.add(Ptemp.transpose()).scalarMultiply(0.5).add(currSC.psdCorr), 1E-6, 1E-16).getL();
			
		while (true)
		{
		    for (int i = 0; i < numStates; i++)
		    {
	    	currSC.sigma.setColumnVector(i, currSC.xhat.add(sqrtP.getColumnVector(i)));
	    	currSC.sigma.setColumnVector(numStates + i, currSC.xhat.subtract(sqrtP.getColumnVector(i)));
		    }
		    
		    double[][] sigData = currSC.sigma.getData();
		    for (int j = 6; j < odCfg.parameters.size() + 6; j++)
		    {
			Settings.Parameter tempep = odCfg.parameters.get(j - 6);
			for (int i = 0; i < numSigmas; i++)
			    sigData[j][i] = FastMath.min(FastMath.max(sigData[j][i], tempep.min), tempep.max);
		    }
		    currSC.sigma.setSubMatrix(sigData, 0, 0);

		    double step = propFinal - propStart;
		    if (odCfg.propStep != 0.0 && (propStart >= bound0 || propFinal >= bound0) &&
			(propStart <= bound1 || propFinal <= bound1))
		    {
			if (odCfg.propStep > 0.0)
			    step = FastMath.min(step, odCfg.propStep);
			else
			    step = FastMath.max(step, odCfg.propStep);
		    }
		    stepSum += step;

		    
		    final double propEnd = propStart + step;
		    
		    
		    if (FastMath.abs(step) > 1.0E-6)
			propagator.propagate(propStart, currSC.sigma, propEnd, currSC.propSigma, false);
		    else
	    	currSC.propSigma.setSubMatrix(currSC.sigma.getData(), 0, 0);

		    
		    currSC.xhatPrev = addColumns(currSC.propSigma).mapMultiplyToSelf(weight);
		    currSC.xhat = new ArrayRealVector(currSC.xhatPrev);
		    currSC.Pprop = odCfg.getProcessNoiseMatrix(stepSum);
		    	
		    

		    
		    
		    for (int i = 0; i < numSigmas; i++)
		    {
			RealVector y = currSC.propSigma.getColumnVector(i).subtract(currSC.xhatPrev);
			currSC.Pprop = currSC.Pprop.add(y.outerProduct(y).scalarMultiply(weight));
			
		    }


		    if (measIndex == odObs.rawMeas.length)
		    {
	    	Estimation.EstimationOutput odout = new Estimation.EstimationOutput();
			odout.time = DataManager.getUTCString(new AbsoluteDate(epoch, propEnd));
			odout.estimatedState = currSC.xhatPrev.toArray();
			odout.propagatedCovariance = currSC.Pprop.getData();
			estOutput.add(odout);
		    }

		    propStart += step;
		    if (FastMath.abs(step) < 1.0E-6 || FastMath.abs(propFinal - propStart) < 1.0E-6)
			break;
		}
		
	    //JPDA here

		
		}

		
		if (measIndex == odObs.rawMeas.length)
		    break EstimationLoop;

		// Lots of JPDA here
		
		
		//update
		for(int objNum = 0; objNum < totalObjNum; objNum++)
		{
			
		SingleObject currSC = Objects.get(objNum);
			
		RealVector rawMeas = null; // same as raw in my code
		for (int i = 0; i < numSigmas; i++)
		{
		    double[] pv = currSC.propSigma.getColumn(i);
		    ssta[0] = new SpacecraftState(new CartesianOrbit(new PVCoordinates(new Vector3D(pv[0], pv[1], pv[2]),
										       new Vector3D(pv[3], pv[4], pv[5])),
								     odCfg.propFrame, tm, Constants.EGM96_EARTH_MU),
						  propagator.getAttitude(tm, pv), odCfg.rsoMass);

		    if (combinedMeas || measNames.length == 1)
		    {
			double[] fitv = odObs.measObjs.get(measIndex).estimate(0, 0, ssta).getEstimatedValue();
			currSC.estimMeas.setColumn(i, fitv);
			if (rawMeas == null)
			    rawMeas = new ArrayRealVector(odObs.measObjs.get(measIndex).getObservedValue());
		    }
		    else
		    {
			double[] fitv = odObs.measObjs.get(measIndex*2).estimate(0, 0, ssta).getEstimatedValue();
			currSC.estimMeas.setEntry(0, i, fitv[0]);
			fitv = odObs.measObjs.get(measIndex*2 + 1).estimate(0, 0, ssta).getEstimatedValue();
			currSC.estimMeas.setEntry(1, i, fitv[0]);
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
			    rawMeas = rawMeas.subtract(new ArrayRealVector(new double[]{currSC.xhatPrev.getEntry(pos),
			    						currSC.xhatPrev.getEntry(biasPos.get(name))}));
			}
			else
			    rawMeas = rawMeas.subtract(new ArrayRealVector(new double[]{currSC.xhatPrev.getEntry(pos)}));
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
		    currSC.xhat = new ArrayRealVector(currSC.xhatPrev.add(odCfg.parameterMatrix.multiply(K).operate(rawMeas.subtract(yhatpre))));
		    currSC.P = currSC.Pprop.subtract(odCfg.parameterMatrix.multiply(K.multiply(Pyy.multiply(K.transpose()))));


		    //For Estimated Parameters, return xhat values to max/min bounds. If far outside the bounds, can result in
		    //Singular Matrix when taking the inverse for smoothing.
		    
		    for (int j = 6; j < odCfg.parameters.size() + 6; j++)
		    {
		    	Settings.Parameter tempep = odCfg.parameters.get(j - 6);
		    	currSC.xhat.setEntry(j, FastMath.min(FastMath.max(currSC.xhat.getEntry(j), tempep.min), tempep.max));
		    }
		    
		    double[] pv = currSC.xhat.toArray();
		    ssta[0] = new SpacecraftState(new CartesianOrbit(new PVCoordinates(new Vector3D(pv[0], pv[1], pv[2]),
										       new Vector3D(pv[3], pv[4], pv[5])),
								     odCfg.propFrame, tm, Constants.EGM96_EARTH_MU),
						  propagator.getAttitude(tm, pv), odCfg.rsoMass);

		    Estimation.EstimationOutput odout = new Estimation.EstimationOutput();
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

		    odout.time = odObs.rawMeas[measIndex].time;
		    odout.station = odObs.rawMeas[measIndex].station;
		    odout.estimatedState = pv;
		    odout.propagatedCovariance = currSC.Pprop.getData();
		    odout.innovationCovariance = Pyy.getData();
		    odout.estimatedCovariance = currSC.P.getData();
		    estOutput.add(odout);
		    
		    
		    //generate post sigma points
		    final RealMatrix Ptemp = currSC.P.scalarMultiply(numStates);
			final RealMatrix sqrtP = new CholeskyDecomposition(
			    Ptemp.add(Ptemp.transpose()).scalarMultiply(0.5).add(currSC.psdCorr), 1E-6, 1E-16).getL();

			Array2DRowRealMatrix postSigma = new Array2DRowRealMatrix(numStates, numSigmas);
			
		    for (int i = 0; i < numStates; i++)
		    {
	    	postSigma.setColumnVector(i, currSC.xhat.add(sqrtP.getColumnVector(i)));
	    	postSigma.setColumnVector(numStates + i, currSC.xhat.subtract(sqrtP.getColumnVector(i)));
		    }

		    double[][] sigData = postSigma.getData();
		    for (int j = 6; j < odCfg.parameters.size() + 6; j++)
		    {
			Settings.Parameter tempep = odCfg.parameters.get(j - 6);
			for (int i = 0; i < numSigmas; i++)
			    sigData[j][i] = FastMath.min(FastMath.max(sigData[j][i], tempep.min), tempep.max);
		    }
		    postSigma.setSubMatrix(sigData, 0, 0);
		    
		    
		    // Store Smoother Data

		    
		    smootherTimeStep smout = new smootherTimeStep();
			
			smout.xpre = MatrixUtils.createColumnRealMatrix(currSC.xhatPrev.toArray());
			smout.xpost = MatrixUtils.createColumnRealMatrix(currSC.xhat.toArray());
			
			smout.Ppre = currSC.Pprop;
			smout.Ppost = currSC.P;
	
			smout.sigPre = MatrixUtils.createRealMatrix(currSC.propSigma.getData());
			smout.sigPost = MatrixUtils.createRealMatrix(postSigma.getData());
	
			smout.tmSmoother = tm;
    		
    		
			if(combinedMeas || measNames.length == 1)
			{
				smout.measObjsSmoother = odObs.measObjs.get(measIndex);
				
			}
			else
			{
				smout.measObjsSmoother = odObs.measObjs.get(2 * measIndex);
				smout.measObjsSmootherNoComb = odObs.measObjs.get(2* measIndex + 1);
			}

			currSC.smootherData.add(smout);
		    
    	}
	    }
		
		
		//Smooth Data
		for(int objNum = 0; objNum < totalObjNum; objNum++)
		{
			
			SingleObject currSC = Objects.get(objNum);
			
			currSC.McReynoldsConsistencyPass = true;
			
			int smSize = currSC.smootherData.size()-1;

			currSC.smootherData.get(smSize).xstar = currSC.smootherData.get(smSize).xpost;
			currSC.smootherData.get(smSize).Pstar = currSC.smootherData.get(smSize).Ppost;			
			

		    for(int i = 0; i < smSize; i++)
		    {
		    	smootherTimeStep smDatak1 = currSC.smootherData.get(smSize - i);
		    	smootherTimeStep smDatak = currSC.smootherData.get(smSize - i - 1);
	
				RealMatrix Csmoother = new Array2DRowRealMatrix(numStates,numStates);
				RealMatrix Asmoother = new Array2DRowRealMatrix(numStates,numStates);
		    	
				for(int j = 0; j < numSigmas; j++) 
				{
					Csmoother = Csmoother.add(((smDatak.sigPost.getColumnMatrix(j).subtract(smDatak.xpost)).multiply( 
					    (smDatak1.sigPre.getColumnMatrix(j).subtract(smDatak1.xpre)).transpose())).scalarMultiply(weight));
				}
				

				Asmoother = Csmoother.multiply(MatrixUtils.inverse(smDatak1.Ppre));
	
				smDatak.xstar = smDatak.xpost.add(Asmoother.multiply(smDatak1.xstar.subtract(smDatak1.xpre)));
				smDatak.Pstar = smDatak.Ppost.add(Asmoother.multiply((smDatak1.Pstar.subtract(smDatak1.Ppre)).multiply(Asmoother.transpose())));
				

				currSC.smootherData.set(smSize - i - 1,smDatak);
				

		    }
		    

		// compute McReynolds Consistency or merge all hypotheses into one hypothesis
	
		double McReynoldsVal = -99999;
		
	    for(int i = 0; i < currSC.smootherData.size()-1;i++)
	    {
	    	smootherTimeStep smDatak1 = currSC.smootherData.get(smSize - i);
	    	smootherTimeStep smDatak = currSC.smootherData.get(smSize - i - 1);

			RealMatrix delx = smDatak.xpost.subtract(smDatak.xstar);	
			RealMatrix delP = smDatak.Ppost.subtract(smDatak.Pstar);
			
			for(int j = 0; j < 5; j++)
			{			
				McReynoldsVal = Math.max(McReynoldsVal, Math.abs(delx.getEntry(j,0)) / Math.sqrt(Math.abs(delP.getEntry(j,j))));
				
				if(Math.abs(delx.getEntry(j,0)) / Math.sqrt(Math.abs(delP.getEntry(j,j))) >= 3)
					currSC.McReynoldsConsistencyPass = false;
			}
		}
		
		
		//store in results
	    
	    SpacecraftState[] smSsta = new SpacecraftState[1];		
	    
	    for(int j = 0; j < currSC.smootherData.size();j++)
	    {
	    	
	    	double[] pv = currSC.smootherData.get(j).xstar.getColumn(0);
		    
		    tm = currSC.smootherData.get(j).tmSmoother;
		    
		    smSsta[0] = new SpacecraftState(new CartesianOrbit(new PVCoordinates(new Vector3D(pv[0], pv[1], pv[2]),
										       new Vector3D(pv[3], pv[4], pv[5])), odCfg.propFrame, tm, Constants.EGM96_EARTH_MU));


	    	//compute smoothed residuals
		    
		    if (combinedMeas || measNames.length == 1)
			{
			    for (int i = 0; i < measNames.length; i++)
			    {
				double[] fitv = currSC.smootherData.get(j).measObjsSmoother.estimate(0, 0, smSsta).getEstimatedValue();
				if (measNames.length == 1)
				{
					estOutput.get(j).postFit.put(measNames[i], fitv);
				}
				else
				{
					estOutput.get(j).postFit.put(measNames[i], new double[] {fitv[i]});

				}
			    }
			}
			else
			{
			    double[] fitv = currSC.smootherData.get(j).measObjsSmoother.estimate(0, 0, smSsta).getEstimatedValue();
			    estOutput.get(j).postFit.put(measNames[0], fitv);
			    if (Rsize > 1)
			    {
				fitv = currSC.smootherData.get(j).measObjsSmootherNoComb.estimate(0, 0, smSsta).getEstimatedValue();
				estOutput.get(j).postFit.put(measNames[1], fitv);
			    }
			}

	    	//store
		    estOutput.get(j).estimatedState = currSC.smootherData.get(j).xstar.getColumn(0);
		    estOutput.get(j).propagatedCovariance = currSC.smootherData.get(j).Pstar.getData();
	    		    	
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
    
	
	
	class smootherTimeStep
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
		ObservedMeasurement measObjsSmootherNoComb;    		
	}
    
	
	
    class SingleObject
    {

	    double[] xInitial;
	    RealMatrix P;
	    final RealMatrix psdCorr;
	    final Array2DRowRealMatrix sigma;
	    final Array2DRowRealMatrix propSigma;
	    final Array2DRowRealMatrix estimMeas;
	    
	    RealMatrix Pprop;
	    ArrayRealVector xhat;
	    RealVector xhatPrev;
	    
    	ArrayList<smootherTimeStep> smootherData;

    	boolean McReynoldsConsistencyPass;
	    
	    public SingleObject(int numStates, int numSigmas, int Rsize)
	    {
	    	
			xInitial = odCfg.getInitialState();
			P = new DiagonalMatrix(odCfg.estmCovariance);
			psdCorr = MatrixUtils.createRealIdentityMatrix(P.getRowDimension()).scalarMultiply(1.0E-6);
			sigma = new Array2DRowRealMatrix(numStates, numSigmas);
			propSigma = new Array2DRowRealMatrix(numStates, numSigmas);
			estimMeas = new Array2DRowRealMatrix(Rsize, numSigmas);
			Pprop = null;
			xhat = new ArrayRealVector(xInitial);
			xhatPrev = new ArrayRealVector(xInitial);
	    	
	    	
	    	Pprop = null;
	    	xhat = new ArrayRealVector(xInitial);
	    	xhatPrev = new ArrayRealVector(xInitial);
	    	smootherData = new ArrayList<smootherTimeStep>();

	    }
	    
		
    }
    }
    
}
