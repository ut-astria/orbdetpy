package org.astria;

import org.hipparchus.geometry.euclidean.threed.Vector3D;
import org.orekit.bodies.GeodeticPoint;
import org.orekit.bodies.OneAxisEllipsoid;
import org.orekit.data.DataSource;
import org.orekit.estimation.measurements.AngularAzEl;
import org.orekit.estimation.measurements.GroundStation;
import org.orekit.estimation.measurements.ObservableSatellite;
import org.orekit.estimation.measurements.modifiers.AngularRadioRefractionModifier;
import org.orekit.files.ccsds.ndm.ParserBuilder;
import org.orekit.files.ccsds.ndm.tdm.ObservationsBlock;
import org.orekit.files.ccsds.ndm.tdm.TdmMetadata;
import org.orekit.files.ccsds.ndm.tdm.TdmParser;
import org.orekit.files.ccsds.section.Segment;
import org.orekit.files.ccsds.utils.FileFormat;
import org.orekit.files.sp3.SP3;
import org.orekit.files.sp3.SP3.SP3Ephemeris;
import org.orekit.files.sp3.SP3Parser;
import org.orekit.frames.FramesFactory;
import org.orekit.frames.Predefined;
import org.orekit.frames.TopocentricFrame;
import org.orekit.gnss.SatelliteSystem;
import org.orekit.models.earth.EarthStandardAtmosphereRefraction;
import org.orekit.propagation.BoundedPropagator;
import org.orekit.propagation.SpacecraftState;
import org.orekit.time.AbsoluteDate;
import org.orekit.time.DateTimeComponents;
import org.orekit.time.GNSSDate;
import org.orekit.time.TimeScalesFactory;
import org.orekit.utils.Constants;
import org.orekit.utils.PVCoordinates;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.nio.file.StandardOpenOption;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;


public final class ReadSP3
{
    public static class SP3OutputMeasJava
    {
        public AbsoluteDate time;
        public String station;
        public double[] values;
        public double[] trueState;
        public SP3OutputMeasJava(AbsoluteDate time, String station, double[] valuesInput, double[] trueStateInput)
        {
            this.time = time;
            this.station = station;
            values = new double[2];
            trueState = new double[6];

            for (int i=0; i<values.length; i++) {
                values[i] = valuesInput[i];
            }

            for (int i=0; i<trueState.length; i++) {
                trueState[i] = trueStateInput[i];
            }

        }

    }
    public static ArrayList<SP3OutputMeasJava> ReadSP3(String tdmFileName, String sp3FileName, String outFilePath, Boolean turnOnLtcorr)
    {

        File filename = new File(tdmFileName);
        String basetdmfilename = filename.getName().split("\\.")[0];
        outFilePath = outFilePath + basetdmfilename + ".json";

        ArrayList<SP3OutputMeasJava> sp3corr_list = new ArrayList<SP3OutputMeasJava>();

        final ObservableSatellite obssat = new ObservableSatellite(0);
        final double[] zeros = new double[]{0.0, 0.0};
        final double[] ones = new double[]{1.0, 1.0};
        final double[] six_ones = new double[]{1.0, 1.0,1.0, 1.0,1.0, 1.0};

        // NORAD ID to PRN data.
        HashMap<String, String> NORAD2PRN = new HashMap<String, String>();


        NORAD2PRN.put("18109A","G04");
        NORAD2PRN.put("00040A","G28");
        NORAD2PRN.put("03005A","G16");


        //Station lookup, station coords in degrees
        HashMap<String, double[]> StationLookup = new HashMap<String, double[]>();

        //Comments show groundstation coordinate reference, then reference frame for measurements.
        StationLookup.put("UTA-ASTRIANet-02",new double[] {32.90305556, 254.4704444, 2225.04}); //WGS, GCRF


        AbsoluteDate referenceDate = null;

        File dir = new File ("C:\\Users\\test\\IdeaProjects\\Validate_GPS_orbits\\orbdetpy\\corrections_test_case\\ASTRIA_DATASET\\test_sp3_corr_feature_input_testing"); //not required now
        File[] directoryListing = dir.listFiles();
        Arrays.sort(directoryListing);

        //for (File child : directoryListing) { - old - used to process for all TDMs in a folder

        String TDMFile =  tdmFileName;
        ArrayList<ArrayList<Measurements.Measurement>> TDMData = Utilities.importTDM(TDMFile, FileFormat.KVN);  // hard-coded value


        TdmParser parser = new ParserBuilder().buildTdmParser();
        Segment<TdmMetadata, ObservationsBlock> temp = parser.parseMessage(new DataSource(TDMFile)).getSegments().get(0);


        // Extract Norad ID and convert to PRN
        String SatNorad  = temp.getMetadata().getParticipants().get(1); // 1 gives participant 1, 2 - participant 2 - station (UTA)
        String StationName = temp.getMetadata().getParticipants().get(2);


        String PRN = NORAD2PRN.get(SatNorad);

        if(PRN == null)
        {
            System.out.println("Unknown Participant");
            //continue;
        }

        // Extract Groundstation and pull lat/long/height
        double[] lla = StationLookup.get(StationName);

        GroundStation gst = null;

        gst = new GroundStation(new TopocentricFrame(new OneAxisEllipsoid(Constants.WGS84_EARTH_EQUATORIAL_RADIUS,
                Constants.WGS84_EARTH_FLATTENING, FramesFactory.getFrame(Predefined.ITRF_CIO_CONV_2010_ACCURATE_EOP)), //DataManager.getFrame("ITRF")
                new GeodeticPoint(lla[0]*Math.PI/180, lla[1]*Math.PI/180, lla[2]), StationName));

// The angular coordinates will be normalized so that the latitude is between ±π/2 and the longitude is between ±π.

        gst.getPrimeMeridianOffsetDriver().setReferenceDate(AbsoluteDate.J2000_EPOCH);
        gst.getPolarOffsetXDriver().setReferenceDate(AbsoluteDate.J2000_EPOCH);
        gst.getPolarOffsetYDriver().setReferenceDate(AbsoluteDate.J2000_EPOCH);


        if (referenceDate == null)
            referenceDate = new AbsoluteDate(DateTimeComponents.parseDateTime(TDMData.get(0).get(0).time.toString()), TimeScalesFactory.getUTC());

        // For each measurement time:
        for(int i = 0; i < TDMData.size(); i++)     // # TDMs
        {

            for(int j = 0; j < TDMData.get(i).size(); j++) {
                // Extract measurement/time

                double RA = TDMData.get(i).get(j).values[0];
                double Dec = TDMData.get(i).get(j).values[1];
                AbsoluteDate tMeas = TDMData.get(i).get(j).time;

                // convert UTC to gps week
                GNSSDate tMeasGPS = new GNSSDate(tMeas, SatelliteSystem.valueOf("GPS"));

                int week = tMeasGPS.getWeekNumber();
                int dayOfWeek = (int) tMeasGPS.getMilliInWeek() / (1000 * 60 * 60 * 24);

                // Read in GPS file corresponding to GPS week

                String SP3File = sp3FileName;

                SP3 SP3ParsedData = null;

                SP3ParsedData = new SP3Parser().parse(new DataSource(SP3File));

                // Pull Ephemeris for correct PRN
                SP3Ephemeris SatelliteEphem = SP3ParsedData.getSatellites().get(PRN);

                if (tMeas.compareTo(SatelliteEphem.getStop()) == 1) {
                    System.out.println(week);
                    System.out.println(dayOfWeek);
                    System.out.println("Skipped Measurement!!");
                    continue;
                }

                // Propagate to measurement time step
                BoundedPropagator SP3Orbit = SatelliteEphem.getPropagator();

                //Output for Orbdetpy ICs
                if (i == 0 && j == 0) {
                    System.out.println(SP3Orbit.getInitialState().getA());
                    System.out.println(SP3Orbit.getInitialState().getE());

                    PVCoordinates PosVel = SP3Orbit.propagate(tMeas).getPVCoordinates(FramesFactory.getEME2000());

                    System.out.println(tMeas.toString() + "," + PosVel.getPosition().getX() + "," + PosVel.getPosition().getY() + "," + PosVel.getPosition().getZ()
                            + "," + PosVel.getVelocity().getX() + "," + PosVel.getVelocity().getY() + "," + PosVel.getVelocity().getZ() + "\n");
                }

                PVCoordinates PosVel = SP3Orbit.propagate(tMeas).getPVCoordinates(FramesFactory.getEME2000());


                SpacecraftState[] sta = new SpacecraftState[]{SP3Orbit.propagate(tMeas)}; //sta[0].getFrame()-GCRF


                // Produce measurement

                double[] obs = null;
                double[] obsAzEl = {0, 0};
                double delRA = 0;
                double delDec = 0;
                double LTcorr = 0;  


                //tMeas = new AbsoluteDate(tMeas, 0.087); //0.1 offset decent
                SpacecraftState[] satState = new SpacecraftState[]{SP3Orbit.propagate(tMeas)};

                // NOTE: Significant difference comes in 'm' level for GCRF and EME_2000

            	/*
        		obs = new AngularRaDec(gst, FramesFactory.getEME2000(), tMeas, zeros, zeros, ones,
					   obssat).estimate(0, 0, satState).getEstimatedValue();
        		*/


                AngularAzEl AzElModel = new AngularAzEl(gst, satState[0].getDate(), zeros, zeros, ones, obssat);
                EarthStandardAtmosphereRefraction refractionModel = new EarthStandardAtmosphereRefraction();
                AngularRadioRefractionModifier refractionCorrection = new AngularRadioRefractionModifier(refractionModel);
                //AzElModel.addModifier(refractionCorrection);
                obsAzEl = AzElModel.estimate(0, 0, satState).getEstimatedValue();

                GeodeticPoint gstPoint = gst.getOffsetGeodeticPoint(sta[0].getDate());
                obs = Conversion.convertAzElToRaDec(tMeas, obsAzEl[0], obsAzEl[1], gstPoint.getLatitude(), gstPoint.getLongitude(), gstPoint.getAltitude(), Predefined.EME2000);

                //LT Travel Time
                // NOTE: If LT corr On automatically saves corr. meas. to an output file
                if (turnOnLtcorr) {
                    LTcorr = sta[0].getPVCoordinates().getPosition().subtract(gst.getBaseFrame().getPVCoordinates(tMeas, sta[0].getFrame()).getPosition()).getNorm() / Constants.SPEED_OF_LIGHT;

                } else {
                    LTcorr = 0;

                }

                AbsoluteDate tMeasLT = new AbsoluteDate(tMeas, -LTcorr);
                SpacecraftState[] satStateLT = new SpacecraftState[]{SP3Orbit.propagate(tMeasLT)};
                Vector3D relPosMinusLT = satStateLT[0].getPVCoordinates(FramesFactory.getEME2000()).getPosition().subtract(gst.getBaseFrame().getPVCoordinates(tMeasLT, FramesFactory.getEME2000()).getPosition());

                //Vector3D aberrationVel = gst.getBaseFrame().getPVCoordinates(tMeasLT, FramesFactory.getICRF()).getVelocity().add(gst.getBaseFrame().getPVCoordinates(tMeasLT, FramesFactory.getEME2000()).getVelocity());
                Vector3D aberrationVel = gst.getBaseFrame().getPVCoordinates(tMeas, FramesFactory.getICRF()).getVelocity();
                // NOTE: Why ICRF and not GCRF?
                delRA = relPosMinusLT.getAlpha() - relPosMinusLT.subtract(aberrationVel.scalarMultiply(LTcorr)).getAlpha();
                delDec = relPosMinusLT.getDelta() - relPosMinusLT.subtract(aberrationVel.scalarMultiply(LTcorr)).getDelta();


                //Bias correction - old
                if (PRN.contentEquals("G04")) {
                    delRA = delRA - 0.532999607 / 180 * Math.PI / 3600;
                    delDec = delDec - 0.293070837 / 180 * Math.PI / 3600;
                } else if (PRN.contentEquals("G28")) {
                    delRA = delRA - 1.171494226 / 180 * Math.PI / 3600;
                    delDec = delDec - 0.193119551 / 180 * Math.PI / 3600;
                }

            	/*
        		AngularAzEl AzElModel = new AngularAzEl(gst, sta[0].getDate(), zeros, zeros, ones, obssat);
				obsAzEl = AzElModel.estimate(0, 0, sta).getEstimatedValue();
        		GeodeticPoint gstPoint = gst.getOffsetGeodeticPoint(sta[0].getDate());

        		double[] refractionRADEC = Conversion.convertAzElToRaDec(tMeas.toString(), obsAzEl[0] , obsAzEl[1] + RefractionCorr.getRefraction(obsAzEl[1]), gstPoint.getLatitude(), gstPoint.getLongitude(), gstPoint.getAltitude(), "EME2000");
        		double[] baseRADEC = Conversion.convertAzElToRaDec(tMeas.toString(), obsAzEl[0] , obsAzEl[1], gstPoint.getLatitude(), gstPoint.getLongitude(), gstPoint.getAltitude(), "EME2000");

        		System.out.println(RefractionCorr.getRefraction(obsAzEl[1]) +"   "+  (baseRADEC[0] - refractionRADEC[0]) +"   "+ (baseRADEC[1] - refractionRADEC[1]));

        		delRA = delRA - (baseRADEC[0] - refractionRADEC[0]);
        		delDec = delDec - (baseRADEC[1] - refractionRADEC[1]);
        		*/


                if (RA - obs[0] > Math.PI) {
                    obs[0] = obs[0] + 2 * Math.PI;
                } else if (RA - obs[0] < -Math.PI) {
                    obs[0] = obs[0] - 2 * Math.PI;
                }

                if (Dec - obs[1] > Math.PI) {
                    obs[1] = obs[1] + 2 * Math.PI;
                } else if (Dec - obs[1] < -Math.PI) {
                    obs[1] = obs[1] - 2 * Math.PI;
                }

                double correctedRA = RA + delRA;
                double correctedDec = Dec + delDec;


                double[] corr_vals = new double[] {correctedRA, correctedDec};
                double[] p = PosVel.getPosition().toArray();
                double[] v = PosVel.getVelocity().toArray();
                double[] trueStateInput = new double[6];

                for (int k = 0; k < trueStateInput.length; k++)
                {
                    if (k < 3)
                        trueStateInput[k] = p[k];
                    else if (k < 6)
                        trueStateInput[k] = v[k - 3];

                }
                SP3OutputMeasJava each_arei = new SP3OutputMeasJava(tMeas, StationName, corr_vals, trueStateInput);
                sp3corr_list.add(each_arei);


                if (turnOnLtcorr) {

                    try {

                        String data = "";

                        if (j==0) {
                            Files.deleteIfExists(new File(outFilePath).toPath());
                        }

                        if (!new File(outFilePath).exists()) {
                            data = "[";
                        }
                        else
                            data = ",\n";


                        data = data + "{\"Time\": \"" + tMeas.toString() + "Z\", \"Station\": \"" +
                                StationName + "\", \"Corr_RA\": " +
                                correctedRA + ", \"Corr_Dec\": " + correctedDec + ", \"NORADID\": \"" + SatNorad + "\"" +
                                ",\"True_State\":[" + PosVel.getPosition().getX() + "," + PosVel.getPosition().getY() + "," + PosVel.getPosition().getZ()
                                + "," + PosVel.getVelocity().getX() + "," + PosVel.getVelocity().getY() + "," + PosVel.getVelocity().getZ() + "]\n" + "}";

                        if (TDMData.get(i).size() - 1 == j)
                            data = data + "]";

                        Files.write(Paths.get(outFilePath), data.getBytes(), StandardOpenOption.CREATE, StandardOpenOption.APPEND);


                    } catch (IOException e) {

                        e.printStackTrace();
                    }


                }


            }

        }

        //} - not using for loop at the moment - old

        return sp3corr_list;


    }

}
