/*
 * Conversion.java - Various conversion functions.
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
import java.util.List;
import org.hipparchus.geometry.euclidean.threed.Vector3D;
import org.hipparchus.linear.MatrixUtils;
import org.hipparchus.linear.RealMatrix;
import org.hipparchus.util.FastMath;
import org.hipparchus.util.MathUtils;
import org.orekit.bodies.GeodeticPoint;
import org.orekit.frames.*;
import org.orekit.time.*;
import org.orekit.utils.CartesianDerivativesFilter;
import org.orekit.utils.PVCoordinates;
import org.orekit.time.TimeScale;
import org.orekit.time.TimeScalesFactory;

public final class Conversion
{
    private Conversion()
    {
    }

    public static double[] transformFrame(Predefined srcFrame, AbsoluteDate time, List<Double> pva, Predefined destFrame)
    {
        Transform xfm = FramesFactory.getFrame(srcFrame).getTransformTo(FramesFactory.getFrame(destFrame), time);
        System.out.println("time from conv.java: " + time);
        System.out.println("size of input: " + pva.size());
        if (pva.size() == 3)
            return(xfm.transformPosition(new Vector3D(pva.get(0), pva.get(1), pva.get(2))).toArray());
        else if (pva.size() == 6)
        {
            PVCoordinates toPv = xfm.transformPVCoordinates(new PVCoordinates(new Vector3D(pva.get(0), pva.get(1), pva.get(2)),
                                                                              new Vector3D(pva.get(3), pva.get(4), pva.get(5))));
            double[] p = toPv.getPosition().toArray();
            double[] v = toPv.getVelocity().toArray();
            return(new double[]{p[0], p[1], p[2], v[0], v[1], v[2]});
        }
        else
        {
            PVCoordinates toPv = xfm.transformPVCoordinates(new PVCoordinates(new Vector3D(pva.get(0), pva.get(1), pva.get(2)),
                                                                              new Vector3D(pva.get(3), pva.get(4), pva.get(5)),
                                                                              new Vector3D(pva.get(6), pva.get(7), pva.get(8))));
            double[] p = toPv.getPosition().toArray();
            double[] v = toPv.getVelocity().toArray();
            double[] a = toPv.getAcceleration().toArray();
            return(new double[]{p[0], p[1], p[2], v[0], v[1], v[2], a[0], a[1], a[2]});
        }
    }

    public static double[] convertAzElToRaDec(AbsoluteDate time, double az, double el, double lat,
                                              double lon, double alt, Predefined frame)
    {
        TopocentricFrame fromFrame = new TopocentricFrame(DataManager.earthShape, new GeodeticPoint(lat, lon, alt), "gs");
        Transform xfm = fromFrame.getTransformTo(FramesFactory.getFrame(frame), time);
        Vector3D toVec = xfm.transformVector(new Vector3D(FastMath.cos(el)*FastMath.sin(az), FastMath.cos(el)*FastMath.cos(az), FastMath.sin(el)));
        double[] xyz = toVec.toArray();
        return(new double[]{MathUtils.normalizeAngle(FastMath.atan2(xyz[1], xyz[0]), FastMath.PI),
                            FastMath.atan2(xyz[2], FastMath.sqrt(xyz[0]*xyz[0]+xyz[1]*xyz[1]))});
    }

    public static double[] convertRaDecToAzEl(Predefined frame, AbsoluteDate time, double ra, double dec,
                                              double lat, double lon, double alt)
    {
        TopocentricFrame toframe = new TopocentricFrame(DataManager.earthShape, new GeodeticPoint(lat, lon, alt), "gs");
        Transform xfm = FramesFactory.getFrame(frame).getTransformTo(toframe, time);
        Vector3D toVec = xfm.transformVector(new Vector3D(FastMath.cos(dec)*FastMath.cos(ra), FastMath.cos(dec)*FastMath.sin(ra), FastMath.sin(dec)));
        double[] xyz = toVec.toArray();
        return(new double[]{MathUtils.normalizeAngle(FastMath.atan2(xyz[0], xyz[1]), FastMath.PI),
                            FastMath.atan2(xyz[2], FastMath.sqrt(xyz[0]*xyz[0]+xyz[1]*xyz[1]))});
    }


    //myEdit to transform covariance from given frame to ICRF
    public static double[] transformFrameCov(Predefined srcFrame, AbsoluteDate time, List<Double> cov, Predefined destFrame)
    {
        // About BEGIN
        // Input:  time, covariance matrix, input frame, destination frame and
        // Output: cov. LTR in dest. frame
        // About END


        // Covariance conversion code - BEGIN
//            // Covariance transformation
////              Ref. Link: https://forum.orekit.org/t/frame-conversion-between-itrf-and-eme2000/242/5
////            To project a covariance matrix from ITRF to EME2000 I suggest you use the method getJacobian 2 from the Transform class.
////
////            If you only use the rotation matrix you will miss a part of the Jacobian.
////
////            Suppose you already have:
////
////            covItrf: a double[][] covariance matrix (6x6) in ITRF frame;
////            date: an AbsoluteDate corresponding to your covariance matrix.
////                    Consider the following Java code:
//


        //final AbsoluteDate oemDate = new AbsoluteDate(2022, 4, 26, 1, 0, 00.000, utc);
        final AbsoluteDate oemDate = time;

        // For testing only - true
        boolean testPrintCovMat = false;

        if (testPrintCovMat) {
            TimeScale utc = TimeScalesFactory.getUTC();
            TimeScale tt = TimeScalesFactory.getTT();
            double[][] covJ2000 = new double[6][6];
            covJ2000[0][0] = 0.009837978166783588; covJ2000[1][0] = -5.271125664375173e-07;
            covJ2000[1][1] = 0.009836306265990186;             covJ2000[2][0] = -1.4410874815939373e-07;
            covJ2000[2][1] = -3.0425483707999703e-07;             covJ2000[2][2] = 0.009837376276431487;
            covJ2000[3][0] = 4.239628801274919e-06;             covJ2000[3][1] =  -2.5602724096419632e-08;
            covJ2000[3][2] =  -6.952791739437874e-09;             covJ2000[3][3] = 1.722940613447208e-08;
            covJ2000[4][0] = -2.6320058752491264e-08;             covJ2000[4][1] = 4.276257759221778e-06;
            covJ2000[4][2] = 1.4682679720170184e-08;             covJ2000[4][3] = 6.927386922572052e-10;
            covJ2000[4][4] = 1.6575005307516436e-08;             covJ2000[5][0] = -7.1445297042062004e-09;
            covJ2000[5][1] = 1.4678260481397842e-08;             covJ2000[5][2] = 4.225900320773944e-06;
            covJ2000[5][3] = 1.881996850866843e-10;             covJ2000[5][4] = -3.1050123167551245e-10;
            covJ2000[5][5] = 1.7636544907952418e-08;

            // filling upper triangular matrix
            for (int i=0; i<6; i++) {
                for (int j=5; j>=i; j--) {
                    covJ2000[i][j] = covJ2000[j][i];
                }
            }

            System.out.println("prinitng cov. matrix in EME2000 frame as given");
            for (double[] row : covJ2000) {
                System.out.println(Arrays.toString(row));
            }

            System.out.println("printing input covariance from Conv.java");
            System.out.println("==========================================");
            for (int i=0; i<cov.size(); i++) {
                System.out.println("i: " + i + "\tcov. value: " + cov.get(i));
            }

        }

        // Form covariance matrix from python input
        double[][] cov_mat = new double[6][6];
        for (int i = 0, k = 0; i < 6; i++)
        {
            for (int j = 0; j < i + 1; j++, k++)
            {
                double value = cov.get(k);
                cov_mat[i][j] = value;
                cov_mat[j][i] = value;
            }
        }

        // pJ2000 contains the covariance matrix in J2000
        final RealMatrix pJ2000 = MatrixUtils.createRealMatrix(cov_mat);

        if (testPrintCovMat) {
            System.out.println("printing input covariance built NOW");
            System.out.println("==========================================");
            for (double[] row : cov_mat) {
                System.out.println(Arrays.toString(row));
            }
        }

        //  Frames definition
        Frame src_frame = FramesFactory.getFrame(srcFrame);  //EME2000
        Frame dest_frame = FramesFactory.getFrame(destFrame); //ICRF


        // METHOD 1 - BELOW
        // Jacobian from ITRF to J2000 at date
        final double[][] jacJ2000ToIcrf = new double[6][6];
        src_frame.getTransformTo(dest_frame, oemDate).getJacobian(CartesianDerivativesFilter.USE_PV, jacJ2000ToIcrf);

        // Covariance transformation, using Hipparchus RealMatrix class to perform the multiplication
        final RealMatrix jJ2000ToIcrf = MatrixUtils.createRealMatrix(jacJ2000ToIcrf);

        // pJ2000 contains the covariance matrix in J2000
        final RealMatrix pIcrf = jJ2000ToIcrf.multiply(pJ2000.multiplyTransposed(jJ2000ToIcrf));
        double[][] matPIcrf = pIcrf.getData();

        if (testPrintCovMat) {
            System.out.println("printing cov in ICRF frame");
            for (double[] row : matPIcrf) {
                System.out.println(Arrays.toString(row));
            }
        }

        // Convert covariance matrix to LTR - below fn from Estimation.java
        int m = pIcrf.getRowDimension();
        double[] out_cov = new double[(int)(0.5*m*(m+1))];
        for (int i = 0, k = 0; i < m; i++)
        {
            for (int j = 0; j <= i; j++)
                out_cov[k++] = pIcrf.getEntry(i, j);
        }
        return(out_cov);
    }

    public static double[] transformFrameCovMethod2(Predefined srcFrame, AbsoluteDate time, List<Double> cov, Predefined destFrame)
    {
        // About BEGIN
        // Input:  time, covariance matrix, input frame, destination frame and
        // Output: cov. LTR in dest. frame
        // About END


        // Covariance conversion code - BEGIN
//            // Covariance transformation
////              Ref. Link: https://forum.orekit.org/t/frame-conversion-between-itrf-and-eme2000/242/5
////            To project a covariance matrix from ITRF to EME2000 I suggest you use the method getJacobian 2 from the Transform class.
////
////            If you only use the rotation matrix you will miss a part of the Jacobian.
////
////            Suppose you already have:
////
////            covItrf: a double[][] covariance matrix (6x6) in ITRF frame;
////            date: an AbsoluteDate corresponding to your covariance matrix.
////                    Consider the following Java code:
//


        //final AbsoluteDate oemDate = new AbsoluteDate(2022, 4, 26, 1, 0, 00.000, utc);
        final AbsoluteDate oemDate = time;

        // For testing only - true
        boolean testPrintCovMat = false;

        if (testPrintCovMat) {
            TimeScale utc = TimeScalesFactory.getUTC();
            TimeScale tt = TimeScalesFactory.getTT();
            double[][] covJ2000 = new double[6][6];
            covJ2000[0][0] = 0.009837978166783588; covJ2000[1][0] = -5.271125664375173e-07;
            covJ2000[1][1] = 0.009836306265990186;             covJ2000[2][0] = -1.4410874815939373e-07;
            covJ2000[2][1] = -3.0425483707999703e-07;             covJ2000[2][2] = 0.009837376276431487;
            covJ2000[3][0] = 4.239628801274919e-06;             covJ2000[3][1] =  -2.5602724096419632e-08;
            covJ2000[3][2] =  -6.952791739437874e-09;             covJ2000[3][3] = 1.722940613447208e-08;
            covJ2000[4][0] = -2.6320058752491264e-08;             covJ2000[4][1] = 4.276257759221778e-06;
            covJ2000[4][2] = 1.4682679720170184e-08;             covJ2000[4][3] = 6.927386922572052e-10;
            covJ2000[4][4] = 1.6575005307516436e-08;             covJ2000[5][0] = -7.1445297042062004e-09;
            covJ2000[5][1] = 1.4678260481397842e-08;             covJ2000[5][2] = 4.225900320773944e-06;
            covJ2000[5][3] = 1.881996850866843e-10;             covJ2000[5][4] = -3.1050123167551245e-10;
            covJ2000[5][5] = 1.7636544907952418e-08;

            // filling upper triangular matrix
            for (int i=0; i<6; i++) {
                for (int j=5; j>=i; j--) {
                    covJ2000[i][j] = covJ2000[j][i];
                }
            }

            System.out.println("prinitng cov. matrix in EME2000 frame as given");
            for (double[] row : covJ2000) {
                System.out.println(Arrays.toString(row));
            }

            System.out.println("printing input covariance from Conv.java");
            System.out.println("==========================================");
            for (int i=0; i<cov.size(); i++) {
                System.out.println("i: " + i + "\tcov. value: " + cov.get(i));
            }

        }

        // Form covariance matrix from python input
        double[][] cov_mat = new double[6][6];
        for (int i = 0, k = 0; i < 6; i++)
        {
            for (int j = 0; j < i + 1; j++, k++)
            {
                double value = cov.get(k);
                cov_mat[i][j] = value;
                cov_mat[j][i] = value;
            }
        }

        // pJ2000 contains the covariance matrix in J2000
        final RealMatrix pJ2000 = MatrixUtils.createRealMatrix(cov_mat);

        if (testPrintCovMat) {
            System.out.println("printing input covariance built NOW");
            System.out.println("==========================================");
            for (double[] row : cov_mat) {
                System.out.println(Arrays.toString(row));
            }
        }

        //  Frames definition
        Frame src_frame = FramesFactory.getFrame(srcFrame);  //EME2000
        Frame dest_frame = FramesFactory.getFrame(destFrame); //ICRF


        // METHOD 2 - BELOW
        // Jacobian from ITRF to J2000 at date
        final double[][] jacJ2000ToIcrf = new double[6][6];
        src_frame.getTransformTo(dest_frame, oemDate).getJacobian(CartesianDerivativesFilter.USE_PV, jacJ2000ToIcrf);

        // Covariance transformation, using Hipparchus RealMatrix class to perform the multiplication
        final RealMatrix jJ2000ToIcrf = MatrixUtils.createRealMatrix(jacJ2000ToIcrf);

        // pJ2000 contains the covariance matrix in J2000
        final RealMatrix pIcrf = jJ2000ToIcrf.multiply(pJ2000.multiplyTransposed(jJ2000ToIcrf));
        double[][] matPIcrf = pIcrf.getData();

        if (testPrintCovMat) {
            System.out.println("printing cov in ICRF frame");
            for (double[] row : matPIcrf) {
                System.out.println(Arrays.toString(row));
            }
        }

        // Convert covariance matrix to LTR - below fn from Estimation.java
        int m = pIcrf.getRowDimension();
        double[] out_cov = new double[(int)(0.5*m*(m+1))];
        for (int i = 0, k = 0; i < m; i++)
        {
            for (int j = 0; j <= i; j++)
                out_cov[k++] = pIcrf.getEntry(i, j);
        }
        return(out_cov);
    }
}
