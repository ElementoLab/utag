from utag.segmentation import utag
import inspect
s = inspect.signature(utag)
assert "return_copy" in s.parameters
assert s.parameters["return_copy"].default is True
print("OK: return_copy param with default True")
