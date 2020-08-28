# Generated by the gRPC Python protocol compiler plugin. DO NOT EDIT!
"""Client and server classes corresponding to protobuf-defined services."""
import grpc

import orbdetpy.rpc.messages_pb2 as messages__pb2


class PropagationStub(object):
    """Missing associated documentation comment in .proto file."""

    def __init__(self, channel):
        """Constructor.

        Args:
            channel: A grpc.Channel.
        """
        self.propagate = channel.unary_unary(
                '/Propagation/propagate',
                request_serializer=messages__pb2.SettingsArray.SerializeToString,
                response_deserializer=messages__pb2.Measurement2DArray.FromString,
                )


class PropagationServicer(object):
    """Missing associated documentation comment in .proto file."""

    def propagate(self, request, context):
        """Missing associated documentation comment in .proto file."""
        context.set_code(grpc.StatusCode.UNIMPLEMENTED)
        context.set_details('Method not implemented!')
        raise NotImplementedError('Method not implemented!')


def add_PropagationServicer_to_server(servicer, server):
    rpc_method_handlers = {
            'propagate': grpc.unary_unary_rpc_method_handler(
                    servicer.propagate,
                    request_deserializer=messages__pb2.SettingsArray.FromString,
                    response_serializer=messages__pb2.Measurement2DArray.SerializeToString,
            ),
    }
    generic_handler = grpc.method_handlers_generic_handler(
            'Propagation', rpc_method_handlers)
    server.add_generic_rpc_handlers((generic_handler,))


 # This class is part of an EXPERIMENTAL API.
class Propagation(object):
    """Missing associated documentation comment in .proto file."""

    @staticmethod
    def propagate(request,
            target,
            options=(),
            channel_credentials=None,
            call_credentials=None,
            insecure=False,
            compression=None,
            wait_for_ready=None,
            timeout=None,
            metadata=None):
        return grpc.experimental.unary_unary(request, target, '/Propagation/propagate',
            messages__pb2.SettingsArray.SerializeToString,
            messages__pb2.Measurement2DArray.FromString,
            options, channel_credentials,
            insecure, call_credentials, compression, wait_for_ready, timeout, metadata)
