# -*- coding: utf-8 -*-
# Generated by the protocol buffer compiler.  DO NOT EDIT!
# source: messages.proto
"""Generated protocol buffer code."""
from google.protobuf import descriptor as _descriptor
from google.protobuf import descriptor_pool as _descriptor_pool
from google.protobuf import message as _message
from google.protobuf import reflection as _reflection
from google.protobuf import symbol_database as _symbol_database
# @@protoc_insertion_point(imports)

_sym_db = _symbol_database.Default()




DESCRIPTOR = _descriptor_pool.Default().AddSerializedFile(b'\n\x0emessages.proto\"!\n\x0cIntegerArray\x12\x11\n\x05\x61rray\x18\x01 \x03(\x05\x42\x02\x10\x01\" \n\x0b\x44oubleArray\x12\x11\n\x05\x61rray\x18\x01 \x03(\x01\x42\x02\x10\x01\",\n\rDouble2DArray\x12\x1b\n\x05\x61rray\x18\x01 \x03(\x0b\x32\x0c.DoubleArray\"H\n\tParameter\x12\r\n\x05value\x18\x01 \x01(\x01\x12\x0b\n\x03min\x18\x02 \x01(\x01\x12\x0b\n\x03max\x18\x03 \x01(\x01\x12\x12\n\nestimation\x18\x04 \x01(\x05\")\n\x05\x46\x61\x63\x65t\x12\x12\n\x06normal\x18\x01 \x03(\x01\x42\x02\x10\x01\x12\x0c\n\x04\x61rea\x18\x02 \x01(\x01\"\x7f\n\x08Maneuver\x12\x0c\n\x04time\x18\x01 \x01(\x01\x12\x15\n\rtrigger_event\x18\x02 \x01(\x05\x12\x1a\n\x0etrigger_params\x18\x03 \x03(\x01\x42\x02\x10\x01\x12\x15\n\rmaneuver_type\x18\x04 \x01(\x05\x12\x1b\n\x0fmaneuver_params\x18\x05 \x03(\x01\x42\x02\x10\x01\"\xad\x01\n\x07Station\x12\x10\n\x08latitude\x18\x01 \x01(\x01\x12\x11\n\tlongitude\x18\x02 \x01(\x01\x12\x10\n\x08\x61ltitude\x18\x03 \x01(\x01\x12\x10\n\x04\x62ias\x18\x04 \x03(\x01\x42\x02\x10\x01\x12\x17\n\x0f\x62ias_estimation\x18\x05 \x01(\x05\x12\x13\n\x0b\x66ov_azimuth\x18\x06 \x01(\x01\x12\x15\n\rfov_elevation\x18\x07 \x01(\x01\x12\x14\n\x0c\x66ov_aperture\x18\x08 \x01(\x01\"8\n\x12MeasurementSetting\x12\x0f\n\x07two_way\x18\x01 \x01(\x08\x12\x11\n\x05\x65rror\x18\x02 \x03(\x01\x42\x02\x10\x01\"\xfc\r\n\x08Settings\x12\x10\n\x08rso_mass\x18\x01 \x01(\x01\x12\x10\n\x08rso_area\x18\x02 \x01(\x01\x12\x1a\n\nrso_facets\x18\x03 \x03(\x0b\x32\x06.Facet\x12 \n\x14rso_solar_array_axis\x18\x04 \x03(\x01\x42\x02\x10\x01\x12\x1c\n\x14rso_solar_array_area\x18\x05 \x01(\x01\x12\x1d\n\x15rso_attitude_provider\x18\x06 \x01(\x05\x12\x1d\n\x11rso_spin_velocity\x18\x07 \x03(\x01\x42\x02\x10\x01\x12!\n\x15rso_spin_acceleration\x18\x08 \x03(\x01\x42\x02\x10\x01\x12\x16\n\x0egravity_degree\x18\t \x01(\x05\x12\x15\n\rgravity_order\x18\n \x01(\x05\x12\x1a\n\x12ocean_tides_degree\x18\x0b \x01(\x05\x12\x19\n\x11ocean_tides_order\x18\x0c \x01(\x05\x12\x16\n\x0ethird_body_sun\x18\r \x01(\x08\x12\x17\n\x0fthird_body_moon\x18\x0e \x01(\x08\x12\x17\n\x0fsolid_tides_sun\x18\x0f \x01(\x08\x12\x18\n\x10solid_tides_moon\x18\x10 \x01(\x08\x12\x12\n\ndrag_model\x18\x11 \x01(\x05\x12$\n\x10\x64rag_coefficient\x18\x12 \x01(\x0b\x32\n.Parameter\x12\'\n\x10\x64rag_MSISE_flags\x18\x13 \x03(\x0b\x32\r.IntegerArray\x12\x15\n\rdrag_exp_rho0\x18\x14 \x01(\x01\x12\x13\n\x0b\x64rag_exp_H0\x18\x15 \x01(\x01\x12\x17\n\x0f\x64rag_exp_Hscale\x18\x16 \x01(\x01\x12\x0e\n\x06rp_sun\x18\x17 \x01(\x08\x12\'\n\x13rp_coeff_reflection\x18\x18 \x01(\x0b\x32\n.Parameter\x12\x1b\n\x13rp_coeff_absorption\x18\x19 \x01(\x01\x12\x1c\n\tmaneuvers\x18\x1a \x03(\x0b\x32\t.Maneuver\x12\x12\n\nprop_start\x18\x1b \x01(\x01\x12\x10\n\x08prop_end\x18\x1c \x01(\x01\x12\x11\n\tprop_step\x18\x1d \x01(\x01\x12\x1e\n\x12prop_initial_state\x18\x1e \x03(\x01\x42\x02\x10\x01\x12\x18\n\x10prop_initial_TLE\x18\x1f \x03(\t\x12\x1b\n\x13prop_inertial_frame\x18  \x01(\t\x12\x1b\n\x13integ_min_time_step\x18! \x01(\x01\x12\x1b\n\x13integ_max_time_step\x18\" \x01(\x01\x12\x1b\n\x13integ_abs_tolerance\x18# \x01(\x01\x12\x1b\n\x13integ_rel_tolerance\x18$ \x01(\x01\x12\x18\n\x10sim_measurements\x18% \x01(\x08\x12)\n\x08stations\x18& \x03(\x0b\x32\x17.Settings.StationsEntry\x12\x31\n\x0cmeasurements\x18\' \x03(\x0b\x32\x1b.Settings.MeasurementsEntry\x12\x1c\n\x10geo_zone_lat_lon\x18( \x03(\x01\x42\x02\x10\x01\x12\x13\n\x0b\x65stm_filter\x18) \x01(\x05\x12\x1b\n\x0f\x65stm_covariance\x18* \x03(\x01\x42\x02\x10\x01\x12\x1e\n\x12\x65stm_process_noise\x18+ \x03(\x01\x42\x02\x10\x01\x12\x1a\n\x12\x65stm_DMC_corr_time\x18, \x01(\x01\x12\x1b\n\x13\x65stm_DMC_sigma_pert\x18- \x01(\x01\x12)\n\x15\x65stm_DMC_acceleration\x18. \x01(\x0b\x32\n.Parameter\x12\x1a\n\x12\x65stm_outlier_sigma\x18/ \x01(\x01\x12\x1b\n\x13\x65stm_outlier_warmup\x18\x30 \x01(\x05\x12 \n\x18\x65stm_smoother_iterations\x18\x31 \x01(\x05\x12\"\n\x1a\x65stm_detection_probability\x18\x32 \x01(\x01\x12\x1f\n\x17\x65stm_gating_probability\x18\x33 \x01(\x01\x12\x1d\n\x15\x65stm_gating_threshold\x18\x34 \x01(\x01\x12\x14\n\x0coutput_flags\x18\x35 \x01(\x05\x12\x12\n\nhyp_sigma1\x18\x36 \x01(\x01\x12\x12\n\nhyp_sigma2\x18\x37 \x01(\x01\x12\x18\n\x10hyp_grid_spacing\x18\x38 \x01(\x01\x12\x13\n\x0bhyp_sma_min\x18\x39 \x01(\x01\x12\x13\n\x0bhyp_sma_max\x18: \x01(\x01\x12\x13\n\x0bhyp_ecc_max\x18; \x01(\x01\x1a\x39\n\rStationsEntry\x12\x0b\n\x03key\x18\x01 \x01(\t\x12\x17\n\x05value\x18\x02 \x01(\x0b\x32\x08.Station:\x02\x38\x01\x1aH\n\x11MeasurementsEntry\x12\x0b\n\x03key\x18\x01 \x01(\x05\x12\"\n\x05value\x18\x02 \x01(\x0b\x32\x13.MeasurementSetting:\x02\x38\x01\")\n\rSettingsArray\x12\x18\n\x05\x61rray\x18\x01 \x03(\x0b\x32\t.Settings\"q\n\x0bMeasurement\x12\x0c\n\x04time\x18\x01 \x01(\x01\x12\x0f\n\x07station\x18\x02 \x01(\t\x12\x12\n\x06values\x18\x03 \x03(\x01\x42\x02\x10\x01\x12\x17\n\x0b\x61ngle_rates\x18\x04 \x03(\x01\x42\x02\x10\x01\x12\x16\n\ntrue_state\x18\x05 \x03(\x01\x42\x02\x10\x01\"/\n\x10MeasurementArray\x12\x1b\n\x05\x61rray\x18\x01 \x03(\x0b\x32\x0c.Measurement\"6\n\x12Measurement2DArray\x12 \n\x05\x61rray\x18\x01 \x03(\x0b\x32\x11.MeasurementArray\"T\n\x13\x44\x65termineOrbitInput\x12\x19\n\x06\x63onfig\x18\x01 \x01(\x0b\x32\t.Settings\x12\"\n\x0cmeasurements\x18\x02 \x03(\x0b\x32\x0c.Measurement\"\xfe\x01\n\x10\x45stimationOutput\x12\x0c\n\x04time\x18\x01 \x01(\x01\x12\x0f\n\x07station\x18\x02 \x01(\t\x12\x1b\n\x0f\x65stimated_state\x18\x03 \x03(\x01\x42\x02\x10\x01\x12!\n\x15propagated_covariance\x18\x04 \x03(\x01\x42\x02\x10\x01\x12!\n\x15innovation_covariance\x18\x05 \x03(\x01\x42\x02\x10\x01\x12 \n\x14\x65stimated_covariance\x18\x06 \x03(\x01\x42\x02\x10\x01\x12\x13\n\x07pre_fit\x18\x07 \x03(\x01\x42\x02\x10\x01\x12\x14\n\x08post_fit\x18\x08 \x03(\x01\x42\x02\x10\x01\x12\x1b\n\x13\x63lutter_probability\x18\t \x01(\x01\"9\n\x15\x45stimationOutputArray\x12 \n\x05\x61rray\x18\x01 \x03(\x0b\x32\x11.EstimationOutput\"V\n\x10MultiTargetInput\x12\x19\n\x06\x63onfig\x18\x01 \x03(\x0b\x32\t.Settings\x12\'\n\x0cmeasurements\x18\x02 \x03(\x0b\x32\x11.MeasurementArray\"\x84\x01\n\x11MultiTargetOutput\x12*\n\nest_output\x18\x01 \x03(\x0b\x32\x16.EstimationOutputArray\x12%\n\x0e\x61ssociated_obs\x18\x02 \x03(\x0b\x32\r.IntegerArray\x12\x1c\n\x10unassociated_obs\x18\x03 \x03(\x05\x42\x02\x10\x01\"{\n\x13TransformFrameInput\x12\x11\n\tsrc_frame\x18\x01 \x01(\t\x12\x10\n\x04time\x18\x02 \x03(\x01\x42\x02\x10\x01\x12\x19\n\x03pva\x18\x03 \x03(\x0b\x32\x0c.DoubleArray\x12\x12\n\ndest_frame\x18\x04 \x01(\t\x12\x10\n\x08UTC_time\x18\x05 \x03(\t\"8\n\x0eImportTDMInput\x12\x11\n\tfile_name\x18\x01 \x01(\t\x12\x13\n\x0b\x66ile_format\x18\x02 \x01(\x05\"\x8d\x01\n\x0b\x41nglesInput\x12\r\n\x05\x66rame\x18\x01 \x01(\t\x12\x10\n\x04time\x18\x02 \x03(\x01\x42\x02\x10\x01\x12\x10\n\x08latitude\x18\x03 \x01(\x01\x12\x11\n\tlongitude\x18\x04 \x01(\x01\x12\x10\n\x08\x61ltitude\x18\x05 \x01(\x01\x12\x12\n\x06\x61ngle1\x18\x06 \x03(\x01\x42\x02\x10\x01\x12\x12\n\x06\x61ngle2\x18\x07 \x03(\x01\x42\x02\x10\x01\"\xa1\x01\n\x19InterpolateEphemerisInput\x12\x14\n\x0csource_frame\x18\x01 \x01(\t\x12\x10\n\x04time\x18\x02 \x03(\x01\x42\x02\x10\x01\x12\x1b\n\x05\x65phem\x18\x03 \x03(\x0b\x32\x0c.DoubleArray\x12\x12\n\nnum_points\x18\x04 \x01(\x05\x12\x12\n\ndest_frame\x18\x05 \x01(\t\x12\x17\n\x0binterp_time\x18\x06 \x03(\x01\x42\x02\x10\x01\x42\x1c\n\x0eorg.astria.rpcB\x08MessagesP\x00\x62\x06proto3')



_INTEGERARRAY = DESCRIPTOR.message_types_by_name['IntegerArray']
_DOUBLEARRAY = DESCRIPTOR.message_types_by_name['DoubleArray']
_DOUBLE2DARRAY = DESCRIPTOR.message_types_by_name['Double2DArray']
_PARAMETER = DESCRIPTOR.message_types_by_name['Parameter']
_FACET = DESCRIPTOR.message_types_by_name['Facet']
_MANEUVER = DESCRIPTOR.message_types_by_name['Maneuver']
_STATION = DESCRIPTOR.message_types_by_name['Station']
_MEASUREMENTSETTING = DESCRIPTOR.message_types_by_name['MeasurementSetting']
_SETTINGS = DESCRIPTOR.message_types_by_name['Settings']
_SETTINGS_STATIONSENTRY = _SETTINGS.nested_types_by_name['StationsEntry']
_SETTINGS_MEASUREMENTSENTRY = _SETTINGS.nested_types_by_name['MeasurementsEntry']
_SETTINGSARRAY = DESCRIPTOR.message_types_by_name['SettingsArray']
_MEASUREMENT = DESCRIPTOR.message_types_by_name['Measurement']
_MEASUREMENTARRAY = DESCRIPTOR.message_types_by_name['MeasurementArray']
_MEASUREMENT2DARRAY = DESCRIPTOR.message_types_by_name['Measurement2DArray']
_DETERMINEORBITINPUT = DESCRIPTOR.message_types_by_name['DetermineOrbitInput']
_ESTIMATIONOUTPUT = DESCRIPTOR.message_types_by_name['EstimationOutput']
_ESTIMATIONOUTPUTARRAY = DESCRIPTOR.message_types_by_name['EstimationOutputArray']
_MULTITARGETINPUT = DESCRIPTOR.message_types_by_name['MultiTargetInput']
_MULTITARGETOUTPUT = DESCRIPTOR.message_types_by_name['MultiTargetOutput']
_TRANSFORMFRAMEINPUT = DESCRIPTOR.message_types_by_name['TransformFrameInput']
_IMPORTTDMINPUT = DESCRIPTOR.message_types_by_name['ImportTDMInput']
_ANGLESINPUT = DESCRIPTOR.message_types_by_name['AnglesInput']
_INTERPOLATEEPHEMERISINPUT = DESCRIPTOR.message_types_by_name['InterpolateEphemerisInput']
IntegerArray = _reflection.GeneratedProtocolMessageType('IntegerArray', (_message.Message,), {
  'DESCRIPTOR' : _INTEGERARRAY,
  '__module__' : 'messages_pb2'
  # @@protoc_insertion_point(class_scope:IntegerArray)
  })
_sym_db.RegisterMessage(IntegerArray)

DoubleArray = _reflection.GeneratedProtocolMessageType('DoubleArray', (_message.Message,), {
  'DESCRIPTOR' : _DOUBLEARRAY,
  '__module__' : 'messages_pb2'
  # @@protoc_insertion_point(class_scope:DoubleArray)
  })
_sym_db.RegisterMessage(DoubleArray)

Double2DArray = _reflection.GeneratedProtocolMessageType('Double2DArray', (_message.Message,), {
  'DESCRIPTOR' : _DOUBLE2DARRAY,
  '__module__' : 'messages_pb2'
  # @@protoc_insertion_point(class_scope:Double2DArray)
  })
_sym_db.RegisterMessage(Double2DArray)

Parameter = _reflection.GeneratedProtocolMessageType('Parameter', (_message.Message,), {
  'DESCRIPTOR' : _PARAMETER,
  '__module__' : 'messages_pb2'
  # @@protoc_insertion_point(class_scope:Parameter)
  })
_sym_db.RegisterMessage(Parameter)

Facet = _reflection.GeneratedProtocolMessageType('Facet', (_message.Message,), {
  'DESCRIPTOR' : _FACET,
  '__module__' : 'messages_pb2'
  # @@protoc_insertion_point(class_scope:Facet)
  })
_sym_db.RegisterMessage(Facet)

Maneuver = _reflection.GeneratedProtocolMessageType('Maneuver', (_message.Message,), {
  'DESCRIPTOR' : _MANEUVER,
  '__module__' : 'messages_pb2'
  # @@protoc_insertion_point(class_scope:Maneuver)
  })
_sym_db.RegisterMessage(Maneuver)

Station = _reflection.GeneratedProtocolMessageType('Station', (_message.Message,), {
  'DESCRIPTOR' : _STATION,
  '__module__' : 'messages_pb2'
  # @@protoc_insertion_point(class_scope:Station)
  })
_sym_db.RegisterMessage(Station)

MeasurementSetting = _reflection.GeneratedProtocolMessageType('MeasurementSetting', (_message.Message,), {
  'DESCRIPTOR' : _MEASUREMENTSETTING,
  '__module__' : 'messages_pb2'
  # @@protoc_insertion_point(class_scope:MeasurementSetting)
  })
_sym_db.RegisterMessage(MeasurementSetting)

Settings = _reflection.GeneratedProtocolMessageType('Settings', (_message.Message,), {

  'StationsEntry' : _reflection.GeneratedProtocolMessageType('StationsEntry', (_message.Message,), {
    'DESCRIPTOR' : _SETTINGS_STATIONSENTRY,
    '__module__' : 'messages_pb2'
    # @@protoc_insertion_point(class_scope:Settings.StationsEntry)
    })
  ,

  'MeasurementsEntry' : _reflection.GeneratedProtocolMessageType('MeasurementsEntry', (_message.Message,), {
    'DESCRIPTOR' : _SETTINGS_MEASUREMENTSENTRY,
    '__module__' : 'messages_pb2'
    # @@protoc_insertion_point(class_scope:Settings.MeasurementsEntry)
    })
  ,
  'DESCRIPTOR' : _SETTINGS,
  '__module__' : 'messages_pb2'
  # @@protoc_insertion_point(class_scope:Settings)
  })
_sym_db.RegisterMessage(Settings)
_sym_db.RegisterMessage(Settings.StationsEntry)
_sym_db.RegisterMessage(Settings.MeasurementsEntry)

SettingsArray = _reflection.GeneratedProtocolMessageType('SettingsArray', (_message.Message,), {
  'DESCRIPTOR' : _SETTINGSARRAY,
  '__module__' : 'messages_pb2'
  # @@protoc_insertion_point(class_scope:SettingsArray)
  })
_sym_db.RegisterMessage(SettingsArray)

Measurement = _reflection.GeneratedProtocolMessageType('Measurement', (_message.Message,), {
  'DESCRIPTOR' : _MEASUREMENT,
  '__module__' : 'messages_pb2'
  # @@protoc_insertion_point(class_scope:Measurement)
  })
_sym_db.RegisterMessage(Measurement)

MeasurementArray = _reflection.GeneratedProtocolMessageType('MeasurementArray', (_message.Message,), {
  'DESCRIPTOR' : _MEASUREMENTARRAY,
  '__module__' : 'messages_pb2'
  # @@protoc_insertion_point(class_scope:MeasurementArray)
  })
_sym_db.RegisterMessage(MeasurementArray)

Measurement2DArray = _reflection.GeneratedProtocolMessageType('Measurement2DArray', (_message.Message,), {
  'DESCRIPTOR' : _MEASUREMENT2DARRAY,
  '__module__' : 'messages_pb2'
  # @@protoc_insertion_point(class_scope:Measurement2DArray)
  })
_sym_db.RegisterMessage(Measurement2DArray)

DetermineOrbitInput = _reflection.GeneratedProtocolMessageType('DetermineOrbitInput', (_message.Message,), {
  'DESCRIPTOR' : _DETERMINEORBITINPUT,
  '__module__' : 'messages_pb2'
  # @@protoc_insertion_point(class_scope:DetermineOrbitInput)
  })
_sym_db.RegisterMessage(DetermineOrbitInput)

EstimationOutput = _reflection.GeneratedProtocolMessageType('EstimationOutput', (_message.Message,), {
  'DESCRIPTOR' : _ESTIMATIONOUTPUT,
  '__module__' : 'messages_pb2'
  # @@protoc_insertion_point(class_scope:EstimationOutput)
  })
_sym_db.RegisterMessage(EstimationOutput)

EstimationOutputArray = _reflection.GeneratedProtocolMessageType('EstimationOutputArray', (_message.Message,), {
  'DESCRIPTOR' : _ESTIMATIONOUTPUTARRAY,
  '__module__' : 'messages_pb2'
  # @@protoc_insertion_point(class_scope:EstimationOutputArray)
  })
_sym_db.RegisterMessage(EstimationOutputArray)

MultiTargetInput = _reflection.GeneratedProtocolMessageType('MultiTargetInput', (_message.Message,), {
  'DESCRIPTOR' : _MULTITARGETINPUT,
  '__module__' : 'messages_pb2'
  # @@protoc_insertion_point(class_scope:MultiTargetInput)
  })
_sym_db.RegisterMessage(MultiTargetInput)

MultiTargetOutput = _reflection.GeneratedProtocolMessageType('MultiTargetOutput', (_message.Message,), {
  'DESCRIPTOR' : _MULTITARGETOUTPUT,
  '__module__' : 'messages_pb2'
  # @@protoc_insertion_point(class_scope:MultiTargetOutput)
  })
_sym_db.RegisterMessage(MultiTargetOutput)

TransformFrameInput = _reflection.GeneratedProtocolMessageType('TransformFrameInput', (_message.Message,), {
  'DESCRIPTOR' : _TRANSFORMFRAMEINPUT,
  '__module__' : 'messages_pb2'
  # @@protoc_insertion_point(class_scope:TransformFrameInput)
  })
_sym_db.RegisterMessage(TransformFrameInput)

ImportTDMInput = _reflection.GeneratedProtocolMessageType('ImportTDMInput', (_message.Message,), {
  'DESCRIPTOR' : _IMPORTTDMINPUT,
  '__module__' : 'messages_pb2'
  # @@protoc_insertion_point(class_scope:ImportTDMInput)
  })
_sym_db.RegisterMessage(ImportTDMInput)

AnglesInput = _reflection.GeneratedProtocolMessageType('AnglesInput', (_message.Message,), {
  'DESCRIPTOR' : _ANGLESINPUT,
  '__module__' : 'messages_pb2'
  # @@protoc_insertion_point(class_scope:AnglesInput)
  })
_sym_db.RegisterMessage(AnglesInput)

InterpolateEphemerisInput = _reflection.GeneratedProtocolMessageType('InterpolateEphemerisInput', (_message.Message,), {
  'DESCRIPTOR' : _INTERPOLATEEPHEMERISINPUT,
  '__module__' : 'messages_pb2'
  # @@protoc_insertion_point(class_scope:InterpolateEphemerisInput)
  })
_sym_db.RegisterMessage(InterpolateEphemerisInput)

if _descriptor._USE_C_DESCRIPTORS == False:

  DESCRIPTOR._options = None
  DESCRIPTOR._serialized_options = b'\n\016org.astria.rpcB\010MessagesP\000'
  _INTEGERARRAY.fields_by_name['array']._options = None
  _INTEGERARRAY.fields_by_name['array']._serialized_options = b'\020\001'
  _DOUBLEARRAY.fields_by_name['array']._options = None
  _DOUBLEARRAY.fields_by_name['array']._serialized_options = b'\020\001'
  _FACET.fields_by_name['normal']._options = None
  _FACET.fields_by_name['normal']._serialized_options = b'\020\001'
  _MANEUVER.fields_by_name['trigger_params']._options = None
  _MANEUVER.fields_by_name['trigger_params']._serialized_options = b'\020\001'
  _MANEUVER.fields_by_name['maneuver_params']._options = None
  _MANEUVER.fields_by_name['maneuver_params']._serialized_options = b'\020\001'
  _STATION.fields_by_name['bias']._options = None
  _STATION.fields_by_name['bias']._serialized_options = b'\020\001'
  _MEASUREMENTSETTING.fields_by_name['error']._options = None
  _MEASUREMENTSETTING.fields_by_name['error']._serialized_options = b'\020\001'
  _SETTINGS_STATIONSENTRY._options = None
  _SETTINGS_STATIONSENTRY._serialized_options = b'8\001'
  _SETTINGS_MEASUREMENTSENTRY._options = None
  _SETTINGS_MEASUREMENTSENTRY._serialized_options = b'8\001'
  _SETTINGS.fields_by_name['rso_solar_array_axis']._options = None
  _SETTINGS.fields_by_name['rso_solar_array_axis']._serialized_options = b'\020\001'
  _SETTINGS.fields_by_name['rso_spin_velocity']._options = None
  _SETTINGS.fields_by_name['rso_spin_velocity']._serialized_options = b'\020\001'
  _SETTINGS.fields_by_name['rso_spin_acceleration']._options = None
  _SETTINGS.fields_by_name['rso_spin_acceleration']._serialized_options = b'\020\001'
  _SETTINGS.fields_by_name['prop_initial_state']._options = None
  _SETTINGS.fields_by_name['prop_initial_state']._serialized_options = b'\020\001'
  _SETTINGS.fields_by_name['geo_zone_lat_lon']._options = None
  _SETTINGS.fields_by_name['geo_zone_lat_lon']._serialized_options = b'\020\001'
  _SETTINGS.fields_by_name['estm_covariance']._options = None
  _SETTINGS.fields_by_name['estm_covariance']._serialized_options = b'\020\001'
  _SETTINGS.fields_by_name['estm_process_noise']._options = None
  _SETTINGS.fields_by_name['estm_process_noise']._serialized_options = b'\020\001'
  _MEASUREMENT.fields_by_name['values']._options = None
  _MEASUREMENT.fields_by_name['values']._serialized_options = b'\020\001'
  _MEASUREMENT.fields_by_name['angle_rates']._options = None
  _MEASUREMENT.fields_by_name['angle_rates']._serialized_options = b'\020\001'
  _MEASUREMENT.fields_by_name['true_state']._options = None
  _MEASUREMENT.fields_by_name['true_state']._serialized_options = b'\020\001'
  _ESTIMATIONOUTPUT.fields_by_name['estimated_state']._options = None
  _ESTIMATIONOUTPUT.fields_by_name['estimated_state']._serialized_options = b'\020\001'
  _ESTIMATIONOUTPUT.fields_by_name['propagated_covariance']._options = None
  _ESTIMATIONOUTPUT.fields_by_name['propagated_covariance']._serialized_options = b'\020\001'
  _ESTIMATIONOUTPUT.fields_by_name['innovation_covariance']._options = None
  _ESTIMATIONOUTPUT.fields_by_name['innovation_covariance']._serialized_options = b'\020\001'
  _ESTIMATIONOUTPUT.fields_by_name['estimated_covariance']._options = None
  _ESTIMATIONOUTPUT.fields_by_name['estimated_covariance']._serialized_options = b'\020\001'
  _ESTIMATIONOUTPUT.fields_by_name['pre_fit']._options = None
  _ESTIMATIONOUTPUT.fields_by_name['pre_fit']._serialized_options = b'\020\001'
  _ESTIMATIONOUTPUT.fields_by_name['post_fit']._options = None
  _ESTIMATIONOUTPUT.fields_by_name['post_fit']._serialized_options = b'\020\001'
  _MULTITARGETOUTPUT.fields_by_name['unassociated_obs']._options = None
  _MULTITARGETOUTPUT.fields_by_name['unassociated_obs']._serialized_options = b'\020\001'
  _TRANSFORMFRAMEINPUT.fields_by_name['time']._options = None
  _TRANSFORMFRAMEINPUT.fields_by_name['time']._serialized_options = b'\020\001'
  _ANGLESINPUT.fields_by_name['time']._options = None
  _ANGLESINPUT.fields_by_name['time']._serialized_options = b'\020\001'
  _ANGLESINPUT.fields_by_name['angle1']._options = None
  _ANGLESINPUT.fields_by_name['angle1']._serialized_options = b'\020\001'
  _ANGLESINPUT.fields_by_name['angle2']._options = None
  _ANGLESINPUT.fields_by_name['angle2']._serialized_options = b'\020\001'
  _INTERPOLATEEPHEMERISINPUT.fields_by_name['time']._options = None
  _INTERPOLATEEPHEMERISINPUT.fields_by_name['time']._serialized_options = b'\020\001'
  _INTERPOLATEEPHEMERISINPUT.fields_by_name['interp_time']._options = None
  _INTERPOLATEEPHEMERISINPUT.fields_by_name['interp_time']._serialized_options = b'\020\001'
  _INTEGERARRAY._serialized_start=18
  _INTEGERARRAY._serialized_end=51
  _DOUBLEARRAY._serialized_start=53
  _DOUBLEARRAY._serialized_end=85
  _DOUBLE2DARRAY._serialized_start=87
  _DOUBLE2DARRAY._serialized_end=131
  _PARAMETER._serialized_start=133
  _PARAMETER._serialized_end=205
  _FACET._serialized_start=207
  _FACET._serialized_end=248
  _MANEUVER._serialized_start=250
  _MANEUVER._serialized_end=377
  _STATION._serialized_start=380
  _STATION._serialized_end=553
  _MEASUREMENTSETTING._serialized_start=555
  _MEASUREMENTSETTING._serialized_end=611
  _SETTINGS._serialized_start=614
  _SETTINGS._serialized_end=2402
  _SETTINGS_STATIONSENTRY._serialized_start=2271
  _SETTINGS_STATIONSENTRY._serialized_end=2328
  _SETTINGS_MEASUREMENTSENTRY._serialized_start=2330
  _SETTINGS_MEASUREMENTSENTRY._serialized_end=2402
  _SETTINGSARRAY._serialized_start=2404
  _SETTINGSARRAY._serialized_end=2445
  _MEASUREMENT._serialized_start=2447
  _MEASUREMENT._serialized_end=2560
  _MEASUREMENTARRAY._serialized_start=2562
  _MEASUREMENTARRAY._serialized_end=2609
  _MEASUREMENT2DARRAY._serialized_start=2611
  _MEASUREMENT2DARRAY._serialized_end=2665
  _DETERMINEORBITINPUT._serialized_start=2667
  _DETERMINEORBITINPUT._serialized_end=2751
  _ESTIMATIONOUTPUT._serialized_start=2754
  _ESTIMATIONOUTPUT._serialized_end=3008
  _ESTIMATIONOUTPUTARRAY._serialized_start=3010
  _ESTIMATIONOUTPUTARRAY._serialized_end=3067
  _MULTITARGETINPUT._serialized_start=3069
  _MULTITARGETINPUT._serialized_end=3155
  _MULTITARGETOUTPUT._serialized_start=3158
  _MULTITARGETOUTPUT._serialized_end=3290
  _TRANSFORMFRAMEINPUT._serialized_start=3292
  _TRANSFORMFRAMEINPUT._serialized_end=3415
  _IMPORTTDMINPUT._serialized_start=3417
  _IMPORTTDMINPUT._serialized_end=3473
  _ANGLESINPUT._serialized_start=3476
  _ANGLESINPUT._serialized_end=3617
  _INTERPOLATEEPHEMERISINPUT._serialized_start=3620
  _INTERPOLATEEPHEMERISINPUT._serialized_end=3781
# @@protoc_insertion_point(module_scope)
