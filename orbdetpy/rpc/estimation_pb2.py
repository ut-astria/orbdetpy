# -*- coding: utf-8 -*-
# Generated by the protocol buffer compiler.  DO NOT EDIT!
# source: estimation.proto
"""Generated protocol buffer code."""
from google.protobuf.internal import builder as _builder
from google.protobuf import descriptor as _descriptor
from google.protobuf import descriptor_pool as _descriptor_pool
from google.protobuf import symbol_database as _symbol_database
# @@protoc_insertion_point(imports)

_sym_db = _symbol_database.Default()


import orbdetpy.rpc.messages_pb2 as messages__pb2


DESCRIPTOR = _descriptor_pool.Default().AddSerializedFile(b'\n\x10\x65stimation.proto\x1a\x0emessages.proto2\xb4\x01\n\nEstimation\x12@\n\x0e\x64\x65termineOrbit\x12\x14.DetermineOrbitInput\x1a\x16.EstimationOutputArray\"\x00\x12\x38\n\rmultiTargetOD\x12\x11.MultiTargetInput\x1a\x12.MultiTargetOutput\"\x00\x12*\n\niodLaplace\x12\x0c.AnglesInput\x1a\x0c.DoubleArray\"\x00\x42%\n\x0eorg.astria.rpcB\x11\x45stimationRequestP\x00\x62\x06proto3')

_builder.BuildMessageAndEnumDescriptors(DESCRIPTOR, globals())
_builder.BuildTopDescriptorsAndMessages(DESCRIPTOR, 'estimation_pb2', globals())
if _descriptor._USE_C_DESCRIPTORS == False:

  DESCRIPTOR._options = None
  DESCRIPTOR._serialized_options = b'\n\016org.astria.rpcB\021EstimationRequestP\000'
  _ESTIMATION._serialized_start=37
  _ESTIMATION._serialized_end=217
# @@protoc_insertion_point(module_scope)
