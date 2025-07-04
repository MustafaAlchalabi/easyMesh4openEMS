default_auto_mesh_hint_for_primitives = {
    'metal_edge_res': None,
    'dirs': 'xyz',
}
default_auto_mesh_hint_for_ports = {
    'is_port': True,
    'metal_edge_res': None,
    'dirs': 'xyz'
 }

def decorate_original_method(method, mesh_dict):
    def decorated_method(*args, **kwargs):
        primitive = method(*args, **kwargs)
        if primitive:
            mesh_dict[primitive] = default_auto_mesh_hint_for_primitives
        return primitive
    return decorated_method

class Material_Copy:
    def __init__(self, original_material, primitives_mesh_setup):
        self.original_material = original_material
        self._mesh_dict = primitives_mesh_setup

    def AddLinPoly(self, *args, **kwargs):#
        decorated_func = decorate_original_method(self.original_material.AddLinPoly, self._mesh_dict)
        return decorated_func(*args, **kwargs)

    def AddBox(self, *args, **kwargs):
        decorated_func = decorate_original_method(self.original_material.AddBox, self._mesh_dict)
        return decorated_func(*args, **kwargs)

    def AddPolygon(self, *args, **kwargs):
        decorated_func = decorate_original_method(self.original_material.AddPolygon, self._mesh_dict)
        return decorated_func(*args, **kwargs)

    def AddCylinder(self, *args, **kwargs):
        decorated_func = decorate_original_method(self.original_material.AddCylinder, self._mesh_dict)
        return decorated_func(*args, **kwargs)

    def AddSphere(self, *args, **kwargs):
        decorated_func = decorate_original_method(self.original_material.AddSphere, self._mesh_dict)
        return decorated_func(*args, **kwargs)

    def AddWire(self, *args, **kwargs):
        decorated_func = decorate_original_method(self.original_material.AddWire, self._mesh_dict)
        return decorated_func(*args, **kwargs)

    def AddCurve(self, *args, **kwargs):
        decorated_func = decorate_original_method(self.original_material.AddCurve, self._mesh_dict)
        return decorated_func(*args, **kwargs)
    
    def AddPoint(self, *args, **kwargs):
        decorated_func = decorate_original_method(self.original_material.AddPoint, self._mesh_dict)
        return decorated_func(*args, **kwargs)

    def __getattr__(self, name):
        return getattr(self.original_material, name)

class CSX_Copy:
    def __init__(self, original_csx, primitives_mesh_setup):
        self.original_csx = original_csx
        self._mesh_dict = primitives_mesh_setup

    def AddMaterial(self, *args, **kwargs):
        material = self.original_csx.AddMaterial(*args, **kwargs)
        return Material_Copy(material, self._mesh_dict)
    
    def AddMetal(self, *args, **kwargs):
        metal = self.original_csx.AddMetal(*args, **kwargs)
        return Material_Copy(metal, self._mesh_dict)

    def __getattr__(self, name):
        return getattr(self.original_csx, name)
        
class FDTD_copy:
    def __init__(self, original_FDTD, primitives_mesh_setup):
        self.original_FDTD = original_FDTD
        self._mesh_dict = primitives_mesh_setup

    def AddLumpedPort(self, *args, **kwargs):
        port = self.original_FDTD.AddLumpedPort(*args, **kwargs)
        if port:
            self._mesh_dict[port] = default_auto_mesh_hint_for_ports
        return port
    
    def AddWaveGuidePort(self, *args, **kwargs):
        port = self.original_FDTD.AddWaveGuidePort(*args, **kwargs)
        if port:
            self._mesh_dict[port] = default_auto_mesh_hint_for_ports
        return port
    
    def AddRectWaveGuidePort(self, *args, **kwargs):
        port = self.original_FDTD.AddRectWaveGuidePort(*args, **kwargs)
        if port:
            self._mesh_dict[port] = default_auto_mesh_hint_for_ports
        return port
    
    def AddMSLPort(self, *args, **kwargs):
        port = self.original_FDTD.AddMSLPort(*args, **kwargs)
        if port:
            self._mesh_dict[port] = default_auto_mesh_hint_for_ports
        return port

    def __getattr__(self, name):
        return getattr(self.original_FDTD, name)
    