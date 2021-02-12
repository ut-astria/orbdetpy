/*
 * CAR.java - Functions for creating a constrained admissible region.
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

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.Arrays;
import java.util.ArrayList;
import org.hipparchus.linear.Array2DRowRealMatrix;
import org.hipparchus.linear.RealMatrix;
import org.hipparchus.analysis.solvers.LaguerreSolver;
import org.hipparchus.complex.Complex;
import org.hipparchus.linear.QRDecomposition;
import org.hipparchus.linear.DecompositionSolver;
import org.orekit.utils.TimeStampedPVCoordinates;
import org.orekit.utils.Constants;

import org.hipparchus.optim.nonlinear.vector.leastsquares.MultivariateJacobianFunction;
import org.hipparchus.optim.nonlinear.vector.leastsquares.ParameterValidator;
import org.hipparchus.optim.nonlinear.vector.leastsquares.LeastSquaresProblem;
import org.hipparchus.optim.nonlinear.vector.leastsquares.LeastSquaresBuilder;
import org.hipparchus.optim.nonlinear.vector.leastsquares.LevenbergMarquardtOptimizer;
import org.hipparchus.optim.nonlinear.vector.leastsquares.LeastSquaresOptimizer;
import org.hipparchus.linear.RealVector;
import org.hipparchus.util.Pair;

import java.io.FileWriter;   // Import the FileWriter class
import java.io.File;  // Import the File class
import java.io.IOException;  // Import the IOException class to handle errors

public class CAR
{
	
	ArrayList <CARGaussianElement> gaussians = new ArrayList<CARGaussianElement>();
	
	/*
	 * Optical CAR
	 * meas1 = RA
	 * meas2 = Dec
	 * meas3 = RaRate
	 * meas4 = DecRate
	 * sigma1 = sigmaRange
	 * sigma2 = sigmaRangeRate
	 */
	
	/*
	 * Range CAR
	 * meas1 = RA
	 * meas2 = Dec
	 * meas3 = Range
	 * meas4 = RangeRate
	 * sigma1 = sigmaRaRate
	 * sigma2 = sigmaDecRate
	 */
	
	public CAR(double meas1, double meas2, double meas3, double meas4, TimeStampedPVCoordinates stationCoords, double sigma1, double sigma2, double gridSpacing, double amin, double amax, double emax, boolean combinedMeas)
    {

		
		double mu = Constants.EGM96_EARTH_MU;
		
		RealMatrix q = new Array2DRowRealMatrix(stationCoords.getPosition().toArray());
		RealMatrix q_d = new Array2DRowRealMatrix(stationCoords.getVelocity().toArray());
		
		// Test Case ////////////////////////////////////////////////////
		meas1 = 5 * Math.PI/180;
		meas2 = -45 * Math.PI/180;
		meas3 = 3600*1000;
		meas4 = -0.6*1000;
		
		double Rearth = Constants.EGM96_EARTH_EQUATORIAL_RADIUS;
		
		q = new Array2DRowRealMatrix(new  double[] {Rearth*Math.cos(Math.PI/3) * Math.cos(0), Rearth*Math.cos(Math.PI/3) * Math.sin(0), Rearth*Math.sin(Math.PI/3)});

		RealMatrix wEarth = new Array2DRowRealMatrix(new  double[] {0, 0, Constants.WGS84_EARTH_ANGULAR_VELOCITY});
		q_d = crossProduct(wEarth,q);

		amin = Rearth;
		amax = 2 * Rearth;
		emax = 0.2;
		
		sigma1 = 1e-6;
		sigma2 = 1e-6;
		gridSpacing = 1e-6;
		
		//////////////////////////////////////////////////////////////////
		
		RealMatrix u_p = new Array2DRowRealMatrix(new  double[] {Math.cos(meas1) * Math.cos(meas2) , Math.sin(meas1) * Math.cos(meas2) , Math.sin(meas2)});
		RealMatrix u_a = new Array2DRowRealMatrix(new  double[] {-1 * Math.sin(meas1) * Math.cos(meas2) , Math.cos(meas1) * Math.cos(meas2) , 0});
		RealMatrix u_d = new Array2DRowRealMatrix(new  double[] {-1 * Math.cos(meas1) * Math.sin(meas2) , -1 * Math.sin(meas1) * Math.sin(meas2) , Math.cos(meas2)});
		
		double w0 = Math.pow(q.getFrobeniusNorm(), 2);
		double w1 = 2 * dotProduct(q_d, u_p);
		double w4 = Math.pow(q_d.getFrobeniusNorm(), 2);
		double w5 = 2 * dotProduct(q, u_p);
		
		if(combinedMeas)
		{
	
			double w2 = Math.pow(meas3, 2) * Math.pow(Math.cos(meas2),2) + Math.pow(meas4, 2);
			double w3 = 2 * meas3 * dotProduct(q_d, u_a) + 2 * meas4 * dotProduct(q_d, u_d);

			
			RealMatrix h1 = crossProduct(q, u_p);
			RealMatrix h2 = crossProduct(u_p, u_a.scalarMultiply(meas3).add(u_d.scalarMultiply(meas4)));
			RealMatrix h3 = crossProduct(u_p, q_d).add(crossProduct(q, u_a.scalarMultiply(meas3).add(u_d.scalarMultiply(meas4))));
			RealMatrix h4 = crossProduct(q, q_d);
			
			
			double c0 = Math.pow(h1.getFrobeniusNorm(), 2);
			double c1 = dotProduct(h1.scalarMultiply(2), h2);
			double c2 = dotProduct(h1.scalarMultiply(2), h3);
			double c3 = dotProduct(h1.scalarMultiply(2), h4);
			double c4 = Math.pow(h2.getFrobeniusNorm(), 2);
			double c5 = dotProduct(h2.scalarMultiply(2), h3);
			double c6 = dotProduct(h2.scalarMultiply(2), h4) + Math.pow(h3.getFrobeniusNorm(), 2);
			double c7 = dotProduct(h3.scalarMultiply(2), h4);
			double c8 = Math.pow(h4.getFrobeniusNorm(), 2);
			
			double[] domain = new double[(int) ((2 * amax)/ gridSpacing)+1];
			
			for(int i = 0; i < domain.length; i++)
			{
				domain[i] =  i * gridSpacing;
			}
			
			Double[] aMinL = new Double[domain.length];
			Double[] aMinU = new Double[domain.length];
			Double[] aMaxL = new Double[domain.length];
			Double[] aMaxU = new Double[domain.length];
			Double[] eMaxL = new Double[domain.length];
			Double[] eMaxU = new Double[domain.length];
			
			
			for(int i = 0; i < domain.length; i++)
			{
	
				double r = domain[i];
				double r2 = r*r;
				double r3 = r2*r;
				double r4 = r3*r;
				
				//compute amin and amax bounds
	
				double F = w2 * r2 + w3 * r + w4 - (2 * mu) / Math.sqrt(r2 + w5 * r + w0);
	
				for(int j =0; j < 2; j++)
				{
					double E;
	
					
					if(j==0){
						E = -mu / (2 * amin);
					}
					else {
						E = -mu / (2 * amax);
					}
					
					if((w1*w1/4 - F + 2 * E) >= 0)
					{
						if (j==0){
							aMinL[i] = -w1/2 - Math.sqrt(w1*w1/4 - F + 2 * E);
							aMinU[i] = -w1/2 + Math.sqrt(w1*w1/4 - F + 2 * E);
							
						}
						else {
							aMaxL[i] = -w1/2 - Math.sqrt(w1*w1/4 - F + 2 * E);
							aMaxU[i] = -w1/2 + Math.sqrt(w1*w1/4 - F + 2 * E);
						}
						
					}
					
				}
				//compute emax bounds
	
				double P = c1 * r2 + c2 * r + c3;
				double U = c4 * r4 + c5 * r3 + c6 * r2 + c7 * r + c8;
				
				double a4 = c0;
				double a3 = P + c0 * w1;
				double a2 = U + c0 * F + w1 * P;
				double a1 = F * P + w1 * U;
				double a0 = F * U + mu*mu * (1 - emax*emax);
				
				LaguerreSolver poly = new LaguerreSolver(1e-12,1e-12);
				
				
				Complex[] roots = poly.solveAllComplex(new double[]{a0,a1,a2,a3,a4},0);
				
				Double[] realRoots = new Double[2];
	
				for(int j = 0; j < roots.length; j ++)
				{
	
					if(Math.abs(roots[j].getImaginary()) < 1e-11)
					{
						if(realRoots[0] == null)
						{
							realRoots[0] = roots[j].getReal();
							realRoots[1] = roots[j].getReal();
						}
						else
						{
							if(realRoots[0] < roots[j].getReal())
							{
								realRoots[0] = roots[j].getReal();
							}
							
							if(realRoots[1] > roots[j].getReal())
							{
								realRoots[1] = roots[j].getReal();
							}
						}
						
					}
				}
	
	
	
				
				if(realRoots[0] != null)
				{
	
					eMaxL[i] = realRoots[1];
					eMaxU[i] = realRoots[0];
	
				}
	
			}	
			
			splitCAR(domain, sigma1, sigma2, aMinU, aMinL, aMaxU, aMaxL, eMaxU, eMaxL);

		}
		else
		{
			
			double a11 = Math.pow(meas3,2)*Math.pow(Math.cos(meas2),2);
			double a13 = meas3 * dotProduct(q_d, u_a);
			double a22 = Math.pow(meas3,2);
			double a23 = meas3 * dotProduct(q_d, u_d);
			double a33 = Math.pow(meas4,2) + w1 * meas4 + w4 - 2 * mu / Math.sqrt(a22 + w5 * meas3 + w0);
			
			RealMatrix r = u_p.scalarMultiply(meas3);

			RealMatrix hp = crossProduct(r, q_d.add(u_p.scalarMultiply(meas4)));
			RealMatrix ha = crossProduct(r.scalarMultiply(meas3), u_a);
			RealMatrix hd = crossProduct(r.scalarMultiply(meas3), u_d);
			
			double b11 = Math.pow(ha.getFrobeniusNorm(), 2);
			double b12 = dotProduct(ha, hd);
			double b22 = Math.pow(hd.getFrobeniusNorm(), 2);
			double b13 = dotProduct(hp, ha);
			double b23 = dotProduct(hp, hd);
			double b33 = Math.pow(hp.getFrobeniusNorm(), 2);

			double c40 = a22*b22;
			double c31 = 2*a22*b12;
			double c30 = 2*a22*b23+2*a23*b22;
			double c22 = a11*b22 + a22*b11;
			double c21 = 2*a22*b13 + 2*a13*b22 + 4*a23*b12;
			double c20 = a22*b33 + 4*a23*b23 + a33*b22;
			double c13 = 2*a11*b12;
			double c12 = 2*a11*b23+4*a13*b12 + 2*a23*b11;
			double c11 = 4*a13*b23 + 4*a23*b13 + 2*a33*b12;
			double c10 = 2*a23*b33 + 2*a33*b23;
			double c04 = a11*b11;
			double c03 = 2*a11*b13 + 2*a13*b11;
			double c02 = a11*b33 + 4*a13*b13 + a33*b11;
			double c01 = 2*a13*b33 + 2 * a33 * b13;
			double c00 = a33*b33;
			
		
			
			//Compute computational domain size for no semi-major axis constraint
			double ac0 = Math.abs(a13/a11);
			double a33bar0 = -a33 + Math.pow(a13,2)/a11 + Math.pow(a23,2)/a22;
			double as0 = Math.abs(a33bar0/a11);
			
			double[] domain = new double[(int) ((2 * (ac0 + as0))/ gridSpacing)+1];
			
			for(int i = 0; i < domain.length; i++)
			{
				domain[i] =  i * gridSpacing - (ac0 + as0);
			}
			
			Double[] aMinL = new Double[domain.length];
			Double[] aMinU = new Double[domain.length];
			Double[] aMaxL = new Double[domain.length];
			Double[] aMaxU = new Double[domain.length];
			Double[] eMaxL = new Double[domain.length];
			Double[] eMaxU = new Double[domain.length];
			

				
			//Fill Interpolation Tables
			double[][] InterpolationTableAMin = new double[180][3];
			double[][] InterpolationTableAMax = new double[180][3];
					
			double a33min = a33 + mu / amin;
			double a33max = a33 + mu / amax;

			double a33barmin = -a33min + Math.pow(a13,2)/a11 + Math.pow(a23,2)/a22;
			double a33barmax = -a33max + Math.pow(a13,2)/a11 + Math.pow(a23,2)/a22;
					
			
			System.out.println(Math.pow(meas4,2));
			System.out.println(w1 * meas4);
			System.out.println(w4);
			System.out.println(2 * mu);
			System.out.println(Math.sqrt(a22 + w5 * meas3 + w0));

			
			System.out.println(a33);
			System.out.println(mu / amin);
			System.out.println(mu / amax);

			
			System.out.println(a33min);
			System.out.println(a33max);
			
			System.out.println(a33barmin);
			System.out.println(a33barmax);

			System.exit(0);
			
			
			double ac = -a13/a11;
			double dc = -a23/a22;
			double asmin = Math.abs(a33barmin/a11);
			double asmax = Math.abs(a33barmax/a11);
			double dsmin = Math.abs(a33barmin/a22);
			double dsmax = Math.abs(a33barmax/a22);	
			
			for(int i = 0; i < 180; i++)
			{
				double theta = (180-i)*Math.PI/180;
				InterpolationTableAMin[i][0] = ac + asmin * Math.cos(theta); //Domain
				InterpolationTableAMin[i][1] = dc - dsmin * Math.sin(theta); //Lower Bound
				InterpolationTableAMin[i][2] = dc + dsmin * Math.sin(theta); // Upper Bound
				
				InterpolationTableAMax[i][0] = ac + asmax * Math.cos(theta); //Domain
				InterpolationTableAMax[i][1] = dc - dsmax * Math.sin(theta); //Lower Bound
				InterpolationTableAMax[i][2] = dc + dsmax * Math.sin(theta); // Upper Bound
			}		
					
			for(int i = 0; i < domain.length; i++)
			{
				double RArate = domain[i];
				//compute amin and amax bounds

				double[] aminData = interpolate(InterpolationTableAMin, RArate);
				double[] amaxData = interpolate(InterpolationTableAMax, RArate);

				if(aminData[0] < aminData[1])
				{
					aMinL[i] = aminData[0];
					aMinU[i] = aminData[1];
				}
				else
				{
					aMinL[i] = new Double(0);
					aMinU[i] = new Double(0);
				}
				

				if(amaxData[0] < amaxData[1])
				{
					aMaxL[i] = amaxData[0];
					aMaxU[i] = amaxData[1];
				}
				else
				{
					aMaxL[i] = new Double(0);
					aMaxU[i] = new Double(0);
				}

				
				//compute emax bounds

				double p0 = c04*Math.pow(RArate,4) + c03*Math.pow(RArate,3) + c02*Math.pow(RArate,2) + c01*RArate + c00 + mu*mu*(1-Math.pow(emax,2));
				double p1 = c13*Math.pow(RArate,3) + c12*Math.pow(RArate,2) + c11*RArate + c10;
				double p2 = c22*Math.pow(RArate,2) + c21*RArate + c20;
				double p3 = c31*RArate + c30;
				double p4 = c40;
				
				LaguerreSolver poly = new LaguerreSolver(1e-12,1e-12);
				
				Complex[] roots = poly.solveAllComplex(new double[]{p0,p1,p2,p3,p4},0);
				
				Double[] realRoots = new Double[2];
	
				
				System.out.println("------------------------------");

				for(int j = 0; j < roots.length; j ++)
				{
					System.out.println(roots[j].getReal() + " + i" + roots[j].getImaginary());
					if(Math.abs(roots[j].getImaginary()) < 1e-11)
					{
						if(realRoots[0] == null)
						{
							realRoots[0] = roots[j].getReal();
							realRoots[1] = roots[j].getReal();
						}
						else
						{
							if(realRoots[0] < roots[j].getReal())
							{
								realRoots[0] = roots[j].getReal();
							}
							
							if(realRoots[1] > roots[j].getReal())
							{
								realRoots[1] = roots[j].getReal();
							}
						}
						
					}
				}
	
				if(realRoots[0] != null)
				{
	
					eMaxL[i] = realRoots[1];
					eMaxU[i] = realRoots[0];
	
				}
				
				
			}

			
			//write gaussians to file
			try {
			      FileWriter myWriter = new FileWriter("RangeCAR.txt");
			      for(int i = 0; i < InterpolationTableAMin.length; i++)
			      {
			          myWriter.write(InterpolationTableAMin[i][0] + "," + InterpolationTableAMin[i][1] + "," + InterpolationTableAMin[i][2] + "," + InterpolationTableAMax[i][0] + "," + InterpolationTableAMax[i][1] + "," + InterpolationTableAMax[i][2]+"\n");
			          //myWriter.write(domain[i] + "," + aMinU[i] + "," + aMinL[i] + "," + aMaxU[i]+ "," + aMaxL[i] + "," + eMaxU[i]+ "," + eMaxL[i] + "\n");
			          //myWriter.write(domain[i] + "," + aMinU[i] + "," + aMinL[i] + "," + aMaxU[i]+ "," + aMaxL[i] + "," + "\n");
			      }
			      myWriter.close();
			    } catch (IOException e) {
			      e.printStackTrace();
			    }
			
			splitCAR(domain, sigma1, sigma2, aMinU, aMinL, aMaxU, aMaxL, eMaxU, eMaxL);

		}
		
		
	    

    }
	
	double[] interpolate(double[][] table, double x)
	{
		int i = 0;
		while(i+2 < table.length && x > table[i+1][0])
		{
			i++;
		}
		
		double slopeL = (table[i][1] - table[i+1][1])/(table[i][0] - table[i+1][0]);
		double slopeU = (table[i][2] - table[i+1][2])/(table[i][0] - table[i+1][0]);
		
		double boundL = slopeL * (x - table[i][0]) + table[i][1];
		double boundU = slopeU * (x - table[i][0]) + table[i][2];
		
		double[] bounds = new double[] {boundL, boundU};
		
		
		return bounds;
		
	}
	
	void splitCAR(double[] domain, double sigma1, double sigma2, Double[] aMinU, Double[] aMinL, Double[] aMaxU, Double[] aMaxL, Double[] eMaxU, Double[] eMaxL)
	{
		double[] GMSplitLibrary = new double[1000];
		GMSplitLibrary = readFile("data/input/uniformSigVals.txt");
		
		// find range of truncated CAR

		int CARIndexStart = 0;
		int CARIndexEnd = 0;
		
		for(int i = 0; i < domain.length; i++)
		{
			
			
			if(aMaxU[i] != null && eMaxU[i] != null && aMinL[i] == null)
			{
				CARIndexStart = i;
				break;
			}
			else if((aMaxU[i] != null && eMaxU[i] != null && aMinL[i] != null))
			{
				if(aMaxU[i] > aMinU[i] && eMaxU[i] > aMinU[i] && aMaxL[i] < aMinL[i] && eMaxL[i] < aMinL[i])
				{
					CARIndexStart = i;
					break;
				}
			}
			
		}

		for(int i = CARIndexStart; i < domain.length; i++)
		{
			if(aMaxU[i] == null || eMaxU[i] == null)
			{
				CARIndexEnd = i-1;
				
				break;
			}
		}		

		
		//build truncated CAR
		 
		ArrayList <boundaryCAR> mainCAR = new ArrayList<boundaryCAR>();
		ArrayList <boundaryCAR> upperCAR = new ArrayList<boundaryCAR>();
		ArrayList <boundaryCAR> lowerCAR = new ArrayList<boundaryCAR>();	
		
		for(int i = 0; i <= CARIndexEnd - CARIndexStart; i++)
		{
						
			if(aMaxU[CARIndexStart + i] != null && eMaxU[CARIndexStart + i] != null && aMinL[CARIndexStart + i] != null)
			{

				if(eMaxU[CARIndexStart + i] > aMinU[CARIndexStart + i])
				{

					boundaryCAR boundaryTemp = new boundaryCAR();
					boundaryTemp.domain = domain[CARIndexStart + i];
					boundaryTemp.lowerBound = aMinU[CARIndexStart + i];
					boundaryTemp.upperBound = Math.min(aMaxU[CARIndexStart + i],eMaxU[CARIndexStart + i]);
					
					upperCAR.add(boundaryTemp);
				}

				if(eMaxL[CARIndexStart + i] < aMinL[CARIndexStart + i])
				{

					boundaryCAR boundaryTemp = new boundaryCAR();
					boundaryTemp.domain = domain[CARIndexStart + i];
					boundaryTemp.lowerBound = Math.max(aMaxL[CARIndexStart + i],eMaxL[CARIndexStart + i]);
					boundaryTemp.upperBound = aMinL[CARIndexStart + i];

					lowerCAR.add(boundaryTemp);
				}
			}
			else if((aMaxU[CARIndexStart + i] != null && eMaxU[CARIndexStart + i] != null && aMinU[CARIndexStart + i] == null))
			{

				boundaryCAR boundaryTemp = new boundaryCAR();
				boundaryTemp.domain = domain[CARIndexStart + i];
				boundaryTemp.lowerBound = Math.max(aMaxL[CARIndexStart + i],eMaxL[CARIndexStart + i]);
				boundaryTemp.upperBound = Math.min(aMaxU[CARIndexStart + i],eMaxU[CARIndexStart + i]);

				
				mainCAR.add(boundaryTemp);
			}
		}
		
		
		
		
		int upperIndex = 0;
		int lowerIndex = 0;
		int mainIndex = 0;
		double Area = 0;
		
		
		for(int i = 0; i < CARIndexEnd - CARIndexStart; i++)
		{
			double currentRho= domain[CARIndexStart + i];
			
			double Area1 = 0;
			double Area2 = 0;
			
			//compute area using trapzeoid rule

			// I dont think i need this area calc, since only relative size of weights matters. Test removing the Area value & see if code still works.
			if(!lowerCAR.isEmpty() && !upperCAR.isEmpty() && upperIndex + 1 == upperCAR.size() && lowerIndex + 1 == lowerCAR.size())
			{

				Area1 += (upperCAR.get(upperIndex).upperBound + mainCAR.get(0).upperBound)*(mainCAR.get(0).domain - upperCAR.get(upperIndex).domain)/2;
				Area2 += (lowerCAR.get(lowerIndex).lowerBound + mainCAR.get(0).lowerBound)*(mainCAR.get(0).domain - lowerCAR.get(lowerIndex).domain)/2;
				
				upperIndex++;
				lowerIndex++;
			}
			
			if(!upperCAR.isEmpty() && upperIndex + 1 < upperCAR.size() && upperCAR.get(upperIndex).domain == currentRho )
			{

				Area1 += (upperCAR.get(upperIndex).upperBound + upperCAR.get(upperIndex + 1).upperBound)*(upperCAR.get(upperIndex + 1).domain - upperCAR.get(upperIndex).domain)/2;
				Area2 += (upperCAR.get(upperIndex).lowerBound + upperCAR.get(upperIndex + 1).lowerBound)*(upperCAR.get(upperIndex + 1).domain - upperCAR.get(upperIndex).domain)/2;
				
				upperIndex++;
			}
			
			if(!lowerCAR.isEmpty() && lowerIndex + 1 < lowerCAR.size() && lowerCAR.get(lowerIndex).domain == currentRho )
			{

				Area1 += (lowerCAR.get(lowerIndex).upperBound + lowerCAR.get(lowerIndex + 1).upperBound)*(lowerCAR.get(lowerIndex + 1).domain - lowerCAR.get(lowerIndex).domain)/2;
				Area2 += (lowerCAR.get(lowerIndex).lowerBound + lowerCAR.get(lowerIndex + 1).lowerBound)*(lowerCAR.get(lowerIndex + 1).domain - lowerCAR.get(lowerIndex).domain)/2;
				
				lowerIndex++;
			}
			

			
			if(!mainCAR.isEmpty() && mainCAR.get(mainIndex).domain == currentRho)
			{

				Area1 += (mainCAR.get(mainIndex).upperBound + mainCAR.get(mainIndex + 1).upperBound)*(mainCAR.get(mainIndex + 1).domain - mainCAR.get(mainIndex).domain)/2;
				Area2 += (mainCAR.get(mainIndex).lowerBound + mainCAR.get(mainIndex + 1).lowerBound)*(mainCAR.get(mainIndex + 1).domain - mainCAR.get(mainIndex).domain)/2;
				
				mainIndex++;
			}
			

			

			Area += Area1 - Area2;

    	}

		upperIndex = 0;
		lowerIndex = 0;
		mainIndex = 0;
		
		
		double[] binLocRho = new double[CARIndexEnd - CARIndexStart];
		double[] binSizeRho = new double[CARIndexEnd - CARIndexStart];
		
		
		for(int i = 0; i < CARIndexEnd - CARIndexStart; i++)
		{
			double currentRho= domain[CARIndexStart + i];
			
			binLocRho[i] = domain[CARIndexStart + i] + (domain[CARIndexStart + i + 1] - domain[CARIndexStart + i]) / 2;
			binSizeRho[i] = 0;

			if(!lowerCAR.isEmpty() && !upperCAR.isEmpty() && upperIndex + 1 == upperCAR.size() && lowerIndex + 1 == lowerCAR.size())
			{
				binSizeRho[i] += ((upperCAR.get(upperIndex).upperBound - lowerCAR.get(lowerIndex).lowerBound)+(mainCAR.get(0).upperBound - mainCAR.get(0).lowerBound))/(2 * Area);
				
				upperIndex++;
				lowerIndex++;
			}
			
			if(!upperCAR.isEmpty() && upperIndex + 1 < upperCAR.size() && upperCAR.get(upperIndex).domain == currentRho )
			{
				binSizeRho[i] += ((upperCAR.get(upperIndex).upperBound - upperCAR.get(upperIndex).lowerBound)+(upperCAR.get(upperIndex+1).upperBound - upperCAR.get(upperIndex+1).lowerBound))/(2 * Area);
				upperIndex++;

			}
			
			if(!lowerCAR.isEmpty() && lowerIndex + 1 < lowerCAR.size() && lowerCAR.get(lowerIndex).domain == currentRho )
			{
				binSizeRho[i] += ((lowerCAR.get(lowerIndex).upperBound - lowerCAR.get(lowerIndex).lowerBound)+(lowerCAR.get(lowerIndex+1).upperBound - lowerCAR.get(lowerIndex+1).lowerBound))/(2 * Area);
				lowerIndex++;

			}
			

			
			if(!mainCAR.isEmpty() && mainCAR.get(mainIndex).domain == currentRho)
			{
				binSizeRho[i] += ((mainCAR.get(mainIndex).upperBound - mainCAR.get(mainIndex).lowerBound)+(mainCAR.get(mainIndex+1).upperBound - mainCAR.get(mainIndex+1).lowerBound))/(2 * Area);
				mainIndex++;
			}

    	}
		
		


			
		
	  	int Jp = -1;
		double sigTemp = sigma1 / (domain[CARIndexEnd] - domain[CARIndexStart]); 
		
		for(int i = 0; i < GMSplitLibrary.length; i++)
		{
			if(sigTemp > GMSplitLibrary[i])
			{
				Jp = i+1;
				break;
			}
			else if(i == (GMSplitLibrary.length-1))
			{
				Jp = GMSplitLibrary.length;
			}
			
		}
		
		
		
		double rangeSigVal = GMSplitLibrary[Jp-1];
		

		double[] rangeMean = new double[Jp];
		double rangeSigma = (domain[CARIndexEnd] - domain[CARIndexStart]) * rangeSigVal;
		
		for(int i=0; i < rangeMean.length; i++)
		{
			rangeMean[i] = domain[CARIndexStart] + (domain[CARIndexEnd] - domain[CARIndexStart]) * (i+1) / (Jp + 1);
		}
		
		
		double[][] gaussianContributions = new double[CARIndexEnd - CARIndexStart][rangeMean.length];
		
		for(int i=0; i < gaussianContributions.length; i++)
		{
			for(int j=0; j < gaussianContributions[i].length; j++)
			{
				gaussianContributions[i][j] = Math.exp(-Math.pow(binLocRho[i]-rangeMean[j], 2) / (2 * rangeSigma * rangeSigma)) / (Math.sqrt(2*  rangeSigma * rangeSigma *  Math.PI)); 
			}
		}
				
		MultivariateJacobianFunction model = new MultivariateJacobianFunction() {
		      public Pair<RealVector, RealMatrix> value(final RealVector point) {

		          RealMatrix jacobian = new Array2DRowRealMatrix(gaussianContributions);
		          RealVector value = jacobian.operate(point);
		          return new Pair<RealVector, RealMatrix>(value, jacobian);
		      }
		};
		
		ParameterValidator boundaries = new ParameterValidator() {
		      public RealVector validate(RealVector params) {

		          for(int i = 0; i < params.getDimension(); i++)
		          {
		        	  if(params.getEntry(i) < 0)
		        		  params.setEntry(i,0);
		        	  
		        	  if(params.getEntry(i) > 1)
		        		  params.setEntry(i,1);
		          }
		          return params;
		      }
		};
		
	   double startArray[] = new double[rangeMean.length]; 
	   Arrays.fill(startArray, 0.0);

		
		
	   LeastSquaresProblem problem = new LeastSquaresBuilder().
                start(startArray).
                model(model).
                target(binSizeRho).
                lazyEvaluation(false).
                maxEvaluations(10000).
                maxIterations(10000).
                parameterValidator(boundaries).
                build();
		
		LeastSquaresOptimizer.Optimum optimum = new LevenbergMarquardtOptimizer().optimize(problem);
		double[] weights = optimum.getPoint().toArray();
		////////////////////////////////////////////////////////////////////////////////////////////
		
		double weightSum = 0;
		
		for(int i=0; i < weights.length; i++)
		{
			weightSum += weights[i];
		}
		
		for(int i=0; i < weights.length; i++)
		{
			weights[i] /= weightSum;
		}
		
		ArrayList <Integer> zeroWeights = new ArrayList<Integer>();

		for(int i =0; i < weights.length; i++)
		{
			if(weights[i] > 1 || weights[i] < 0)
			zeroWeights.add(i);
			//throw(new RuntimeException(String.format("Invalid Gaussian weight %f", weights[i])));

			if(weights[i] == 0)
			zeroWeights.add(i);
		}


		// make arraylist for mean, std, weights? and CARGaussianElement Object?
		for(int i = 0; i < rangeMean.length; i++)
		{
			//skip weights that are <= 0 or >1
			if(zeroWeights.contains(i))
			{
				continue;
			}
			ArrayList <CARGaussianElement> tempGaussians = new ArrayList<CARGaussianElement>();

			double currentRho = -1;
			
			for(int j = 0; j < CARIndexEnd - CARIndexStart; j++)
			{
				if(rangeMean[i] >= domain[CARIndexStart + j] && rangeMean[i] <= domain[CARIndexStart + j + 1])
				{
					currentRho = domain[CARIndexStart + j];
					break;

				}
				
			}
			
			
			boolean skipLower = false;
			boolean skipMain = false;

			//upper loop
			if(!upperCAR.isEmpty())
			{
				for(int k = 0; k < upperCAR.size(); k++)
				{
					if(upperCAR.get(k).domain == currentRho)
					{
						
						if(k + 1 < upperCAR.size())
						{
							
							tempGaussians.addAll(generateHypotheses(upperCAR, k, sigma2, rangeMean[i], rangeSigma, GMSplitLibrary));
							
							skipMain = true;
						}
						else
						{
							//construct temp object for the region where change from upper/lower to main occurs
							
							ArrayList <boundaryCAR> tempCAR = new ArrayList<boundaryCAR>();

							tempCAR.add(upperCAR.get(k));
							tempCAR.get(0).lowerBound = lowerCAR.get(lowerCAR.size() - 1).lowerBound;
							tempCAR.add(mainCAR.get(0));
							
							
							
							tempGaussians.addAll(generateHypotheses(tempCAR, 0, sigma2,rangeMean[i], rangeSigma, GMSplitLibrary));

							skipLower = true;
							
						}

					}
					
				}
			}
			
			//lower loop
			if(!lowerCAR.isEmpty() && !skipLower)
			{
				for(int k = 0; k < lowerCAR.size(); k++)
				{
					if(lowerCAR.get(k).domain == currentRho)
					{

						tempGaussians.addAll(generateHypotheses(lowerCAR, k, sigma2,rangeMean[i], rangeSigma, GMSplitLibrary));
						
						skipMain = true;
											
						
					}
					
				}
			}
			
			//main loop
			if(!skipMain)
			{
				for(int k = 0; k < mainCAR.size(); k++)
				{
					if(mainCAR.get(k).domain == currentRho)
					{

						tempGaussians.addAll(generateHypotheses(mainCAR, k, sigma2,rangeMean[i], rangeSigma, GMSplitLibrary));
						
						skipMain = true;
											
						
					}
					
				}
			}
			
			
			for(int k = 0; k < tempGaussians.size(); k++)
			{
				
				tempGaussians.get(k).weight = weights[i] / tempGaussians.size();
			}
			
			
			
			gaussians.addAll(tempGaussians);
    	}
		
		
		//write gaussians to file
		try {
		      FileWriter myWriter = new FileWriter("Gaussians.txt");
		      System.out.println(gaussians.size());
		      for(int i = 0; i < gaussians.size(); i++)
		      {
		          myWriter.write(gaussians.get(i).rangeMean + "," + gaussians.get(i).rangeRateMean + "," + gaussians.get(i).rangeStd + "," + gaussians.get(i).rangeRateStd+ "," + gaussians.get(i).weight+"\n");
		      }
		      myWriter.close();
		    } catch (IOException e) {
		      e.printStackTrace();
		    }
	}
	
	ArrayList <CARGaussianElement> getCAR()
	{
		return gaussians;
	}
	
	RealMatrix crossProduct(RealMatrix vect_A, RealMatrix vect_B) 
	{ 
		RealMatrix cross_P = new Array2DRowRealMatrix(3,1);
		
	    cross_P.setEntry(0,0, vect_A.getEntry(1,0) * vect_B.getEntry(2,0) - vect_A.getEntry(2,0) * vect_B.getEntry(1,0)); 
	    cross_P.setEntry(1,0, vect_A.getEntry(2,0) * vect_B.getEntry(0,0) - vect_A.getEntry(0,0) * vect_B.getEntry(2,0)); 
	    cross_P.setEntry(2,0, vect_A.getEntry(0,0) * vect_B.getEntry(1,0) - vect_A.getEntry(1,0) * vect_B.getEntry(0,0)); 
	    
	    return cross_P;
	}
	
	
	
	double dotProduct(RealMatrix vect_A, RealMatrix vect_B) 
	{ 
		double dot_p = vect_A.getEntry(0,0) * vect_B.getEntry(0,0) + vect_A.getEntry(1,0) * vect_B.getEntry(1,0) + vect_A.getEntry(2,0) * vect_B.getEntry(2,0);
	    
	    return dot_p;
	}

	
	
	
	ArrayList <CARGaussianElement> generateHypotheses(ArrayList <boundaryCAR> currentCAR, int k, double sigma2, double rangeMean, double rangeSigma, double[] GMSplitLibrary)
	{
		ArrayList <CARGaussianElement> tempGaussians = new ArrayList<CARGaussianElement>();

		//interpolate for min/max
		
		double max_rr = (currentCAR.get(k + 1).upperBound - currentCAR.get(k).upperBound) / (currentCAR.get(k + 1).domain - currentCAR.get(k).domain) 
						* (rangeMean - currentCAR.get(k).domain) + currentCAR.get(k).upperBound;
		
		double min_rr = (currentCAR.get(k + 1).lowerBound - currentCAR.get(k).lowerBound) / (currentCAR.get(k + 1).domain - currentCAR.get(k).domain) 
						* (rangeMean - currentCAR.get(k).domain) + currentCAR.get(k).lowerBound;
		
		//find number of mean points (Jp)
				
	  	int Jp_rr = -1;
		double sigTemp = sigma2 / (max_rr - min_rr); 
		
		for(int w = 0; w < GMSplitLibrary.length; w++)
		{
			if(sigTemp > GMSplitLibrary[w])
			{
				Jp_rr = w+1;
				break;
			}
			else if(w == (GMSplitLibrary.length-1))
			{
				Jp_rr = GMSplitLibrary.length;
			}
		}

		//generate mean points
		double rangeRateSigma = (max_rr - min_rr) * GMSplitLibrary[Jp_rr-1];
		

		for(int w=0; w < Jp_rr; w++)
		{
			CARGaussianElement temp = new CARGaussianElement();
			
			
			temp.rangeMean = rangeMean;
			temp.rangeRateMean = min_rr + (max_rr - min_rr) * (w+1) / (Jp_rr + 1);
			temp.rangeRateStd = rangeRateSigma;
			temp.rangeStd = rangeSigma;
			
			tempGaussians.add(temp);
			
		}
		
		return tempGaussians;
		
	}
	
		
	
	double[] readFile(String filename)
	{
		double[] splitLib = new double[1000];


		//read sigma library
		BufferedReader objReader = null;

		try {
	
		int counter = 0;
		String strCurrentLine;
	   
		objReader = new BufferedReader(new FileReader(filename));
	
		while ((strCurrentLine = objReader.readLine()) != null) {   
		splitLib[counter] = Double.valueOf(strCurrentLine);
		counter++;
		}
	
		} catch (IOException e) {
	
		e.printStackTrace();
	
		} finally {
		
			try {
			if (objReader != null)
				objReader.close();
			} catch (IOException ex) {
				ex.printStackTrace();
			}
		}
		  
		  return splitLib;
	}
	
	class boundaryCAR
	{
    	double domain;
    	double upperBound;
    	double lowerBound;
		
	}
	
	
	class CARGaussianElement
	{
		
		double rangeMean;
		double rangeRateMean;
		double rangeStd;
		double rangeRateStd;
		double weight;
		
		
	}
	
}
