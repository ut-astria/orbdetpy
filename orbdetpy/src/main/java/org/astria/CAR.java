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

import java.io.FileWriter;   // Import the FileWriter class
import java.io.File;  // Import the File class
import java.io.IOException;  // Import the IOException class to handle errors

public class CAR
{
	
	ArrayList <CARGaussianElement> gaussians = new ArrayList<CARGaussianElement>();
	
	public CAR(double ra, double ra_d, double dec, double dec_d, TimeStampedPVCoordinates stationCoords, double sigma_range, double sigma_rr, double rangeSample, double amin, double amax, double emax)
    {

		/*
		//this is in rad, also its RA & dec not az el
		double ra = 0.87406;
		double ra_d = 7.29274E-5;
		double dec = -0.07313;
		double dec_d = 6.81829E-7;
		
		RealMatrix q = new Array2DRowRealMatrix(new  double[] {5655719, 1925262, 2234277});
		RealMatrix q_d = new Array2DRowRealMatrix(new  double[] {-140.4001, 412.1821, 0.2260});	
		
		double sigma_range = 15000;
		double sigma_rr = 60;
		
		double amin = 30000000 * 0.9;
		double amax = 47000000 * 1.1;
		double emax = 0.7;
		
		double rangeSample = 5000; 
		*/
		
		double mu = Constants.EGM96_EARTH_MU;

		RealMatrix q = new Array2DRowRealMatrix(stationCoords.getPosition().toArray());
		RealMatrix q_d = new Array2DRowRealMatrix(stationCoords.getVelocity().toArray());	
		
		RealMatrix u_p = new Array2DRowRealMatrix(new  double[] {Math.cos(ra) * Math.cos(dec) , Math.sin(ra) * Math.cos(dec) , Math.sin(dec)});
		RealMatrix u_a = new Array2DRowRealMatrix(new  double[] {-1 * Math.sin(ra) * Math.cos(dec) , Math.cos(ra) * Math.cos(dec) , 0});
		RealMatrix u_d = new Array2DRowRealMatrix(new  double[] {-1 * Math.cos(ra) * Math.sin(dec) , -1 * Math.sin(ra) * Math.sin(dec) , Math.cos(dec)});

		double w0 = Math.pow(q.getFrobeniusNorm(), 2);
		double w1 = 2 * dotProduct(q_d, u_p);
		double w2 = Math.pow(ra_d, 2) * Math.pow(Math.cos(dec),2) + Math.pow(dec_d, 2);
		double w3 = 2 * ra_d * dotProduct(q_d, u_a) + 2 * dec_d * dotProduct(q_d, u_d);
		double w4 = Math.pow(q_d.getFrobeniusNorm(), 2);
		double w5 = 2 * dotProduct(q, u_p);

		RealMatrix h1 = crossProduct(q, u_p);
		RealMatrix h2 = crossProduct(u_p, u_a.scalarMultiply(ra_d).add(u_d.scalarMultiply(dec_d)));
		RealMatrix h3 = crossProduct(u_p, q_d).add(crossProduct(q, u_a.scalarMultiply(ra_d).add(u_d.scalarMultiply(dec_d))));
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
		
		double[] rho = new double[(int) ((2 * amax)/ rangeSample)+1];
		
		for(int i = 0; i < rho.length; i++)
		{
			rho[i] =  i * rangeSample;
		}
		
		Double[] rhoRate_amin_l = new Double[rho.length];
		Double[] rhoRate_amin_u = new Double[rho.length];
		Double[] rhoRate_amax_l = new Double[rho.length];
		Double[] rhoRate_amax_u = new Double[rho.length];
		Double[] rhoRate_emax_l = new Double[rho.length];
		Double[] rhoRate_emax_u = new Double[rho.length];
		
		double[] GMSplitLibrary = new double[1000];
		
		GMSplitLibrary = readFile("data/input/uniformSigVals.txt");
		
		for(int i = 0; i < rho.length; i++)
		{


			double r = rho[i];
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
						rhoRate_amin_l[i] = -w1/2 - Math.sqrt(w1*w1/4 - F + 2 * E);
						rhoRate_amin_u[i] = -w1/2 + Math.sqrt(w1*w1/4 - F + 2 * E);
						
					}
					else {
						rhoRate_amax_l[i] = -w1/2 - Math.sqrt(w1*w1/4 - F + 2 * E);
						rhoRate_amax_u[i] = -w1/2 + Math.sqrt(w1*w1/4 - F + 2 * E);
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


				
				rhoRate_emax_l[i] = realRoots[1];
				rhoRate_emax_u[i] = realRoots[0];
			}

		}		

		// find range of truncated CAR

		int CARIndexStart = 0;
		int CARIndexEnd = 0;
		
		for(int i = 0; i < rho.length; i++)
		{
			
			
			if(rhoRate_amax_u[i] != null && rhoRate_emax_u[i] != null && rhoRate_amin_l[i] == null)
			{
				CARIndexStart = i;
				break;
			}
			else if((rhoRate_amax_u[i] != null && rhoRate_emax_u[i] != null && rhoRate_amin_l[i] != null))
			{
				if(rhoRate_amax_u[i] > rhoRate_amin_u[i] && rhoRate_emax_u[i] > rhoRate_amin_u[i] && rhoRate_amax_l[i] < rhoRate_amin_l[i] && rhoRate_emax_l[i] < rhoRate_amin_l[i])
				{
					CARIndexStart = i;
					break;
				}
			}
			
		}

		for(int i = CARIndexStart; i < rho.length; i++)
		{
			if(rhoRate_amax_u[i] == null || rhoRate_emax_u[i] == null)
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
						

			
			
			if(rhoRate_amax_u[CARIndexStart + i] != null && rhoRate_emax_u[CARIndexStart + i] != null && rhoRate_amin_l[CARIndexStart + i] != null)
			{

				//////////////////////////////////////////////////////////////////////////////////////////// THIS NEEDS TO BE VERIFIED. USED IMAGE ON JAH PAPER PG 3 TO GENERATE LOGIC
				//////////////////////////////////////////////////////////////////////////////////////////// I think it works, checked by using matlab, changed emax to 0.2
				if(rhoRate_emax_u[CARIndexStart + i] > rhoRate_amin_u[CARIndexStart + i])
				{

					boundaryCAR boundaryTemp = new boundaryCAR();
					boundaryTemp.rho = rho[CARIndexStart + i];
					boundaryTemp.lowerBound = rhoRate_amin_u[CARIndexStart + i];
					boundaryTemp.upperBound = Math.min(rhoRate_amax_u[CARIndexStart + i],rhoRate_emax_u[CARIndexStart + i]);
					
					upperCAR.add(boundaryTemp);
				}

				if(rhoRate_emax_l[CARIndexStart + i] < rhoRate_amin_l[CARIndexStart + i])
				{

					boundaryCAR boundaryTemp = new boundaryCAR();
					boundaryTemp.rho = rho[CARIndexStart + i];
					boundaryTemp.lowerBound = Math.max(rhoRate_amax_l[CARIndexStart + i],rhoRate_emax_l[CARIndexStart + i]);
					boundaryTemp.upperBound = rhoRate_amin_l[CARIndexStart + i];

					lowerCAR.add(boundaryTemp);
				}
			}
			else if((rhoRate_amax_u[CARIndexStart + i] != null && rhoRate_emax_u[CARIndexStart + i] != null && rhoRate_amin_u[CARIndexStart + i] == null))
			{

				boundaryCAR boundaryTemp = new boundaryCAR();
				boundaryTemp.rho = rho[CARIndexStart + i];
				boundaryTemp.lowerBound = Math.max(rhoRate_amax_l[CARIndexStart + i],rhoRate_emax_l[CARIndexStart + i]);
				boundaryTemp.upperBound = Math.min(rhoRate_amax_u[CARIndexStart + i],rhoRate_emax_u[CARIndexStart + i]);

				mainCAR.add(boundaryTemp);
			}
		}
		
		
		
		
		int upperIndex = 0;
		int lowerIndex = 0;
		int mainIndex = 0;
		double Area = 0;
		
		
		for(int i = 0; i < CARIndexEnd - CARIndexStart; i++)
		{
			double currentRho= rho[CARIndexStart + i];
			
			double Area1 = 0;
			double Area2 = 0;
			
			//compute area using trapzeoid rule

			// I dont think i need this area calc, since only relative size of weights matters. Test removing the Area value & see if code still works.
			if(!lowerCAR.isEmpty() && !upperCAR.isEmpty() && upperIndex + 1 == upperCAR.size() && lowerIndex + 1 == lowerCAR.size())
			{

				Area1 += (upperCAR.get(upperIndex).upperBound + mainCAR.get(0).upperBound)*(mainCAR.get(0).rho - upperCAR.get(upperIndex).rho)/2;
				Area2 += (lowerCAR.get(lowerIndex).lowerBound + mainCAR.get(0).lowerBound)*(mainCAR.get(0).rho - lowerCAR.get(lowerIndex).rho)/2;
				
				upperIndex++;
				lowerIndex++;
			}
			
			if(!upperCAR.isEmpty() && upperIndex + 1 < upperCAR.size() && upperCAR.get(upperIndex).rho == currentRho )
			{

				Area1 += (upperCAR.get(upperIndex).upperBound + upperCAR.get(upperIndex + 1).upperBound)*(upperCAR.get(upperIndex + 1).rho - upperCAR.get(upperIndex).rho)/2;
				Area2 += (upperCAR.get(upperIndex).lowerBound + upperCAR.get(upperIndex + 1).lowerBound)*(upperCAR.get(upperIndex + 1).rho - upperCAR.get(upperIndex).rho)/2;
				
				upperIndex++;
			}
			
			if(!lowerCAR.isEmpty() && lowerIndex + 1 < lowerCAR.size() && lowerCAR.get(lowerIndex).rho == currentRho )
			{

				Area1 += (lowerCAR.get(lowerIndex).upperBound + lowerCAR.get(lowerIndex + 1).upperBound)*(lowerCAR.get(lowerIndex + 1).rho - lowerCAR.get(lowerIndex).rho)/2;
				Area2 += (lowerCAR.get(lowerIndex).lowerBound + lowerCAR.get(lowerIndex + 1).lowerBound)*(lowerCAR.get(lowerIndex + 1).rho - lowerCAR.get(lowerIndex).rho)/2;
				
				lowerIndex++;
			}
			

			
			if(!mainCAR.isEmpty() && mainCAR.get(mainIndex).rho == currentRho)
			{

				Area1 += (mainCAR.get(mainIndex).upperBound + mainCAR.get(mainIndex + 1).upperBound)*(mainCAR.get(mainIndex + 1).rho - mainCAR.get(mainIndex).rho)/2;
				Area2 += (mainCAR.get(mainIndex).lowerBound + mainCAR.get(mainIndex + 1).lowerBound)*(mainCAR.get(mainIndex + 1).rho - mainCAR.get(mainIndex).rho)/2;
				
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
			double currentRho= rho[CARIndexStart + i];
			
			binLocRho[i] = rho[CARIndexStart + i] + (rho[CARIndexStart + i + 1] - rho[CARIndexStart + i]) / 2;
			binSizeRho[i] = 0;

			if(!lowerCAR.isEmpty() && !upperCAR.isEmpty() && upperIndex + 1 == upperCAR.size() && lowerIndex + 1 == lowerCAR.size())
			{
				binSizeRho[i] += ((upperCAR.get(upperIndex).upperBound - lowerCAR.get(lowerIndex).lowerBound)+(mainCAR.get(0).upperBound - mainCAR.get(0).lowerBound))/(2 * Area);
				
				upperIndex++;
				lowerIndex++;
			}
			
			if(!upperCAR.isEmpty() && upperIndex + 1 < upperCAR.size() && upperCAR.get(upperIndex).rho == currentRho )
			{
				binSizeRho[i] += ((upperCAR.get(upperIndex).upperBound - upperCAR.get(upperIndex).lowerBound)+(upperCAR.get(upperIndex+1).upperBound - upperCAR.get(upperIndex+1).lowerBound))/(2 * Area);
				upperIndex++;

			}
			
			if(!lowerCAR.isEmpty() && lowerIndex + 1 < lowerCAR.size() && lowerCAR.get(lowerIndex).rho == currentRho )
			{
				binSizeRho[i] += ((lowerCAR.get(lowerIndex).upperBound - lowerCAR.get(lowerIndex).lowerBound)+(lowerCAR.get(lowerIndex+1).upperBound - lowerCAR.get(lowerIndex+1).lowerBound))/(2 * Area);
				lowerIndex++;

			}
			

			
			if(!mainCAR.isEmpty() && mainCAR.get(mainIndex).rho == currentRho)
			{
				binSizeRho[i] += ((mainCAR.get(mainIndex).upperBound - mainCAR.get(mainIndex).lowerBound)+(mainCAR.get(mainIndex+1).upperBound - mainCAR.get(mainIndex+1).lowerBound))/(2 * Area);
				mainIndex++;
			}

    	}
		
		


			
		
	  	int Jp = -1;
		double sigTemp = sigma_range / (rho[CARIndexEnd] - rho[CARIndexStart]); 
		
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
		double rangeSigma = (rho[CARIndexEnd] - rho[CARIndexStart]) * rangeSigVal;
		
		for(int i=0; i < rangeMean.length; i++)
		{
			rangeMean[i] = rho[CARIndexStart] + (rho[CARIndexEnd] - rho[CARIndexStart]) * (i+1) / (Jp + 1);
		}
		
		
		double[][] gaussianContributions = new double[CARIndexEnd - CARIndexStart][rangeMean.length];
		
		for(int i=0; i < gaussianContributions.length; i++)
		{
			for(int j=0; j < gaussianContributions[i].length; j++)
			{
				gaussianContributions[i][j] = Math.exp(-Math.pow(binLocRho[i]-rangeMean[j], 2) / (2 * rangeSigma * rangeSigma)) / (Math.sqrt(2*  rangeSigma * rangeSigma *  Math.PI)); // SOME OF THESE VALUES SLIGHTLY OFF FROM MATLAB, NOT SURE THE REASON FOR DIFFERENCE.
			}
		}
		
		double[] weights = constrainedLeastSquares(new Array2DRowRealMatrix(gaussianContributions), new Array2DRowRealMatrix(binSizeRho), 1, 0);

		double weightSum = 0;
		
		for(int i=0; i < weights.length; i++)
		{
			weightSum += weights[i];
		}
		
		for(int i=0; i < weights.length; i++)
		{
			weights[i] /= weightSum;
		}
		
		
		for(int i =0; i < weights.length; i++)
		{
			if(weights[i] > 1 || weights[i] < 0)
			throw(new RuntimeException(String.format("Invalid Gaussian weight %f", weights[i])));

		}


		// make arraylist for mean, std, weights? and CARGaussianElement Object?
		for(int i = 0; i < rangeMean.length; i++)
		{
			ArrayList <CARGaussianElement> tempGaussians = new ArrayList<CARGaussianElement>();

			double currentRho = -1;
			
			for(int j = 0; j < CARIndexEnd - CARIndexStart; j++)
			{
				if(rangeMean[i] >= rho[CARIndexStart + j] && rangeMean[i] <= rho[CARIndexStart + j + 1])
				{
					currentRho = rho[CARIndexStart + j];
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
					if(upperCAR.get(k).rho == currentRho)
					{
						
						if(k + 1 < upperCAR.size())
						{
							
							tempGaussians.addAll(generateHypotheses(upperCAR, k, sigma_rr, rangeMean[i], rangeSigma, GMSplitLibrary));
							
							skipMain = true;
						}
						else
						{
							//construct temp object for the region where change from upper/lower to main occurs
							
							ArrayList <boundaryCAR> tempCAR = new ArrayList<boundaryCAR>();

							tempCAR.add(upperCAR.get(k));
							tempCAR.get(0).lowerBound = lowerCAR.get(lowerCAR.size() - 1).lowerBound;
							tempCAR.add(mainCAR.get(0));
							
							
							
							tempGaussians.addAll(generateHypotheses(tempCAR, 0, sigma_rr,rangeMean[i], rangeSigma, GMSplitLibrary));

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
					if(lowerCAR.get(k).rho == currentRho)
					{

						tempGaussians.addAll(generateHypotheses(lowerCAR, k, sigma_rr,rangeMean[i], rangeSigma, GMSplitLibrary));
						
						skipMain = true;
											
						
					}
					
				}
			}
			
			//main loop
			if(!skipMain)
			{
				for(int k = 0; k < mainCAR.size(); k++)
				{
					if(mainCAR.get(k).rho == currentRho)
					{

						tempGaussians.addAll(generateHypotheses(mainCAR, k, sigma_rr,rangeMean[i], rangeSigma, GMSplitLibrary));
						
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
		
		/*
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
		    */

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

	

	ArrayList <CARGaussianElement> generateHypotheses(ArrayList <boundaryCAR> currentCAR, int k, double sigma_rr, double rangeMean, double rangeSigma, double[] GMSplitLibrary)
	{
		ArrayList <CARGaussianElement> tempGaussians = new ArrayList<CARGaussianElement>();

		//interpolate for min/max
		
		double max_rr = (currentCAR.get(k + 1).upperBound - currentCAR.get(k).upperBound) / (currentCAR.get(k + 1).rho - currentCAR.get(k).rho) 
						* (rangeMean - currentCAR.get(k).rho) + currentCAR.get(k).upperBound;
		
		double min_rr = (currentCAR.get(k + 1).lowerBound - currentCAR.get(k).lowerBound) / (currentCAR.get(k + 1).rho - currentCAR.get(k).rho) 
						* (rangeMean - currentCAR.get(k).rho) + currentCAR.get(k).lowerBound;
		
		//find number of mean points (Jp)
				
	  	int Jp_rr = -1;
		double sigTemp = sigma_rr / (max_rr - min_rr); 
		
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
	
	
	//Solve min(|Ax - b|) for x given bounds on x (upper and lower)
	//Uses the following algorithm https://www.stat.berkeley.edu/~stark/Preprints/bvls.pdf
	//QR decomposition used
	double[] constrainedLeastSquares(RealMatrix A, RealMatrix b, double UB, double LB)
	{
		ArrayList <Integer> free = new ArrayList<Integer>();
		ArrayList <Integer> lower = new ArrayList<Integer>();
		ArrayList <Integer> upper = new ArrayList<Integer>();

		ArrayList <Integer> boundComponents = new ArrayList<Integer>(); // Bound components checks to make sure variable that is just bound is not freed again.

		
		boolean roundError = false;
		int tstar = - 99;
		double tstarVal;
		char tstarOrigBound = 'a';
		
		RealMatrix x = new Array2DRowRealMatrix(A.getColumnDimension(),1);
		
		
		
		for(int  i = 0; i < A.getColumnDimension(); i++)
		{
			lower.add(i);
			x.setEntry(i,0, LB);
		}

		mainLoop:
		while(true)
		{
			System.out.println("Free: " + free.size() + "   Lower: " + lower.size() + "   Upper: " + upper.size());
			// Round Error Test. Component becomes free then moves out of free region.
			if(roundError == true)
			{
				
				System.out.println("Round Error Adjustment");
				
				free.remove(Integer.valueOf(tstar));
				
				if(tstarOrigBound == 'l')
				lower.add(tstar);
				else
				upper.add(tstar);
				
			}
			
			
			RealMatrix w = (A.transpose()).multiply(b.subtract(A.multiply(x)));
			
			//Set problem component term in w to zero
			if(roundError == true)
			{
				roundError = false;
				w.setEntry(tstar,0,0);
			}
			
			//Step 2 test for Kuhn Tucker convergence.
			boolean finished = true;
					
			for(int i = 0; i < lower.size(); i++)
			{
				if(w.getEntry(lower.get(i),0) >= 0)
				{
					finished = false;
				}
			}
			
			for(int i = 0; i < upper.size(); i++)
			{
				if(w.getEntry(upper.get(i),0) >= 0)
				{
					finished = false;
				}
			}
			//
			
			if(b.getFrobeniusNorm() * 1e-12 > b.subtract(A.multiply(x)).getFrobeniusNorm() )
			{
				System.out.println("Solution Essentially Optimal");
				finished = true;
			}
			
			if(free.size() == A.getColumnDimension() || finished)
			{
				break mainLoop;
			}
			else
			{
				//step 4
				/*
				if(lower.size() > 0)
				{
					tstar = lower.get(0);
					tstarVal = w.getEntry(tstar,0);
					tstarOrigBound = 'l';
				}
				else
				{
					tstar = upper.get(0);
					tstarVal = -w.getEntry(tstar,0);
					tstarOrigBound = 'u';

				}*/
				
				tstarVal = Double.NEGATIVE_INFINITY;
				
				//System.out.println(free.size());
				
				for(int i = 0; i < lower.size(); i++)
				{
					if(w.getEntry(lower.get(i),0) > tstarVal && !boundComponents.contains(Integer.valueOf(lower.get(i))))
					{
						tstar = lower.get(i);
						tstarVal = w.getEntry(tstar,0);
						tstarOrigBound = 'l';
					}
				}
				
				for(int i = 0; i < upper.size(); i++)
				{
					if(-w.getEntry(upper.get(i),0) > tstarVal && !boundComponents.contains(Integer.valueOf(upper.get(i))))
					{
						tstar = upper.get(i);
						tstarVal = -w.getEntry(tstar,0);
						tstarOrigBound = 'u';
					}
				}
				
				boundComponents = new ArrayList<Integer>();

				
				//step 5
				lower.remove(Integer.valueOf(tstar));
				upper.remove(Integer.valueOf(tstar));
				
				free.add(tstar);

				//step 6
				RealMatrix bprime = b.copy();
				
				RealMatrix Atemp = A.copy();
				RealMatrix Aprime = new Array2DRowRealMatrix(A.getRowDimension(),free.size());
				
				int[] jprime = new int[free.size()];

				for(int i =0; i < free.size(); i++)
				{
					Atemp.setColumn(free.get(i), new double[A.getRowDimension()]);
					Aprime.setColumn(i, A.getColumn(free.get(i)));
					
					
				}
				
				bprime = bprime.subtract(Atemp.multiply(x));
				
				DecompositionSolver solver = new QRDecomposition(Aprime).getSolver();

				RealMatrix z = solver.solve(bprime);
				
				ArrayList <Integer> zOutOfBounds = new ArrayList<Integer>();

				//step 7
				for(int i = 0; i < z.getRowDimension(); i++)
				{
					if(z.getEntry(i,0) < LB || z.getEntry(i,0) > UB)
					{
						zOutOfBounds.add(i);
					}
				}
				
				if(zOutOfBounds.size() == 0)
				{
					for(int i = 0; i < z.getRowDimension(); i++)
					{
						x.setEntry(free.get(i), 0, z.getEntry(i,0));
					}
					
					continue mainLoop;
				}
				else
				{
					
					//step 8-9.
					
					double alpha = Math.min((LB - x.getEntry(free.get(zOutOfBounds.get(0)),0))/(z.getEntry(zOutOfBounds.get(0),0) - x.getEntry(free.get(zOutOfBounds.get(0)),0)), 
										   (UB - x.getEntry(free.get(zOutOfBounds.get(0)),0))/(z.getEntry(zOutOfBounds.get(0),0) - x.getEntry(free.get(zOutOfBounds.get(0)),0)));
					
					for(int i = 1; i < zOutOfBounds.size(); i++)
					{
						
							double alphaTemp = Math.min((LB - x.getEntry(free.get(zOutOfBounds.get(i)),0))/(z.getEntry(zOutOfBounds.get(i),0) - x.getEntry(free.get(zOutOfBounds.get(i)),0)), 
														(UB - x.getEntry(free.get(zOutOfBounds.get(i)),0))/(z.getEntry(zOutOfBounds.get(i),0) - x.getEntry(free.get(zOutOfBounds.get(i)),0)));
						
							alpha = Math.min(alpha, alphaTemp);
						
					}
					
					for(int i = 0; i < free.size(); i++)
					{		
						x.setEntry(free.get(i),0, x.getEntry(free.get(i),0) + alpha * (z.getEntry(i,0) - x.getEntry(free.get(i),0)));
							
					}
					
					if((x.getEntry(free.size()-1,0) >= UB && tstarOrigBound == 'u') || (x.getEntry(free.size()-1,0) <= LB && tstarOrigBound == 'l'))
					{
						
						//roundError = true;
						//continue mainLoop;
						
					}
					
					for(int i = 0; i < free.size(); i++)
					{		
						if(x.getEntry(free.get(i),0) >= UB)
						{
							int tempIndex = free.get(i);
							
							x.setEntry(tempIndex,0,UB);

							free.remove(Integer.valueOf(tempIndex));
							
							upper.add(tempIndex);
							boundComponents.add(tempIndex);

							
						}
						else if(x.getEntry(free.get(i),0) <= LB)
						{
							int tempIndex = free.get(i);
									
							x.setEntry(tempIndex,0,LB);

							free.remove(Integer.valueOf(tempIndex));
							
							lower.add(tempIndex);
							boundComponents.add(tempIndex);
						}
							
					}
					
					
					
				}
				
			}
	
		}

		
		return x.getColumn(0);

		
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
    	double rho;
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
