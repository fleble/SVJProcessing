from typing import Optional, List, Dict
import numpy as np
import tritonclient.grpc as triton_grpc
import tritonclient.http as triton_http
import torch

# tools required to do inference of pytorch model on GPU using triton


# for converting torch model to proper format needed by triton using jit
def convert_torch_for_triton(torch_model, jit_model_path=None):

  if jit_model_path == None:
    jit_model_path = '/srv/utils/data/GNNTagger/jit_model.pt'

  try:
    jit_model = torch.jit.script(torch_model)
    torch.jit.save(jit_model, jit_model_path)
    print(f"""Torch model was successfully converted for triton inference \
and can be found at {jit_model_path}""")
    return torch.jit.load(jit_model_path)
  except:  #jit became stable in torch >=1.8
    source_command = '\'source /cvmfs/sft.cern.ch/lcg/views/LCG_102/x86_64-centos7-gcc8-opt/setup.sh\''
    raise Exception(f"""TritonUtils error: Torch model had problems with the \
conversion, perhaps your torch version is too old. You can source a newer \
version of, torch for the conversion with {source_command}""")


# for interacting with triton server
# modeled off of https://github.com/rkansal47/HHbbVV/blob/main/src/HHbbVV/postprocessing/triton_check.py
class wrapped_triton:
  def __init__(self, model_url: str, ) -> None:
    fullprotocol, location = model_url.split("://")
    _, protocol = fullprotocol.split("+")
    address, model, version = location.split("/")

    self._protocol = protocol
    self._address = address
    self._model = model
    self._version = version

    # check connection to server, throw error if connection doesn't work
    if self._protocol == "grpc":
      self._client = triton_grpc.InferenceServerClient(url=self._address,
                                                       verbose=False,
                                                       ssl=True)
      self._triton_protocol = triton_grpc
    elif self._protocol == "http":
      self._client = triton_http.InferenceServerClient(url=self._address,
                                                       verbose=False,
                                                       concurrency=12,
                                                       )
      self._triton_protocol = triton_http
    else:
      raise ValueError(
          f"{self._protocol} does not encode a valid protocol (grpc or http)")

  def __call__(self, input_dict: Dict[str, np.ndarray]) -> np.ndarray:
    '''
    Run inference of model on triton server
    '''

    # put inputs in proper format
    inputs = []
    for key in input_dict:
      input = self._triton_protocol.InferInput(key, input_dict[key].shape,
                                               "FP32")
      input.set_data_from_numpy(input_dict[key])
      inputs.append(input)

    output = self._triton_protocol.InferRequestedOutput("softmax__0")

    # make request to server for inference
    request = self._client.infer(self._model,
                                 model_version=self._version,
                                 inputs=inputs,
                                 outputs=[output],
                                 )
    out = request.as_numpy("softmax__0")

    return out