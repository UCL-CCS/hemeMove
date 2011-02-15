import numpy as N

from vtk import vtkPlaneWidget, vtkPolyData, vtkPolyDataMapper, \
     vtkActor

from ..Util.Observer import Observable, ObservableList, NotifyOptions

import pdb
class PlacedIoletList(ObservableList):
    def __init__(self, *args, **kwargs):
        ObservableList.__init__(self, *args, **kwargs)

        # Add observers to add/remove observers for items. Note that
        # for removal we must trigger the action BEFORE the removal
        # happens, so that we can get the item to remove the observer.
        self.AddObserver("@INSERTION", self.HandleInsertion)
        self.AddObserver("@REMOVAL", self.HandlePreRemoval,
                         options=NotifyOptions(BEFORE_CHANGE=True,
                                               AFTER_CHANGE=False))
        return
    
    def SetItemEnabledChangeHandler(self, handler):
        self._ItemEnabledChangeHandler = handler
        return
    
    def SetInteractor(self, iact):
        self.Interactor = iact
        return
    
    def HandleInsertion(self, change):
        self[change.index].AddObserver('Enabled', self._ItemEnabledChangeHandler)
        self[change.index].widget.SetInteractor(self.Interactor)
        return
    
    def HandlePreRemoval(self, change):
        self[change.index].RemoveObserver('Enabled', self._ItemEnabledChangeHandler)
        return
    pass

class PlacedIolet(Observable):
    
    def __init__(self):
        # Working arrays, constructed once for speed.
        self._c = N.zeros(3)
        self._n = N.zeros(3)
        self._p1 = N.zeros(3)
        self._p2 = N.zeros(3)

        self.widget = vtkPlaneWidget()
        self.widget.SetRepresentationToOutline()
        # planeWidget.PlaceWidget(clickPos)
        
        self.representation = vtkPolyData()
        # This is effectively a copy and is guaranteed to be up to
        # date when InteractionEvent or EndInteraction events are
        # invoked
        self.widget.AddObserver("InteractionEvent", self._SyncRepresentation)
        self.widget.GetPolyData(self.representation)

        self.mapper = vtkPolyDataMapper()
        self.mapper.SetInput(self.representation)
        
        self.actor = vtkActor()
        self.actor.SetMapper(self.mapper)
        self.actor.GetProperty().SetColor(self.colour)
        
        self.Enabled = False
        
        self.AddObserver('Enabled', self.EnabledSet)
#        self.AddDependency('IsCentreValid', 'Centre')
#        self.AddDependency('IsNormalValid', 'Normal')
#        self.AddDependency('IsRadiusValid', 'Radius')
#        self.AddDependency('CanShow', 'IsCentreValid')
#        self.AddDependency('CanShow', 'IsNormalValid')
#        self.AddDependency('CanShow', 'IsRadiusValid')
        
    #     self.iolet.AddObserver('Centre.@ANY', self.MyIoletCentreChanged)
    #     self.iolet.AddObserver('Normal.@ANY', self.MyIoletNormalChanged)
    #     self.iolet.AddObserver('Radius', self.MyIoletRadiusChanged)
    #     return
    # def MyIoletCentreChanged(self, change):
    #     self.SetCentre(self.iolet.Centre.x,
    #                    self.iolet.Centre.y,
    #                    self.iolet.Centre.z)
        
    def EnabledSet(self, change):
        if self.widget.GetInteractor() is None:
            return
        
        if self.Enabled:
            self.widget.On()
        else:
            self.widget.Off()
            pass
        return
    
    def _SyncRepresentation(self, obj, evt):
        self.widget.GetPolyData(self.representation)
        return
    
    def SetCentre(self, centre):
        assert N.alltrue(N.isfinite(centre))
        self.widget.SetCenter(centre)
        # Force the plane to be updated
#        self.widget.InvokeEvent("InteractionEvent")
        return
    def GetCentre(self):
        return self.widget.GetCenter()
    Centre = property(GetCentre, SetCentre)
#    @property
#    def IsCentreValid(self):
#        return N.alltrue(N.isfinite(self.Centre))
    
    def SetNormal(self, normal):
        self.widget.SetNormal(normal)
        # Force the plane to be updated
#        self.widget.InvokeEvent("InteractionEvent")
        return
    def GetNormal(self):
        return self.widget.GetNormal()
    Normal = property(GetNormal, SetNormal)
#    @property
#    def IsNormalValid(self):
#        normal = self.Normal
#        return N.alltrue(N.isfinite(normal)) and \
#               N.dot(normal, normal) >1e-6

    _v1 = N.array([0., 0., 1.])
    _v2 = N.array([0., 1., 0.])
    def CalcFirstPlaneUnitVector(self):
        e1 = N.cross(self._n, self._v1)
        e1Sq = N.dot(e1, e1)
        
        e2 = N.cross(self._n, self._v2)
        e2Sq = N.dot(e2, e2)
        
        if e1Sq > e2Sq:
            self._p1 = e1 / N.sqrt(e1Sq)
        else:
            self._p1 = e2 / N.sqrt(e2Sq)
            pass
        return
    
    def CalcSecondPlaneUnitVector(self):
        self._p2 = N.cross(self._n, self._p1)
        return
    
    def SetRadius(self, radius):
        pdb.set_trace()
        # Get into numpy vectors
        self.widget.GetCenter(self._c)
        self.widget.GetNormal(self._n)
        # Calc unit vectors
        self.CalcFirstPlaneUnitVector()
        self.CalcSecondPlaneUnitVector()
        # Scale
        self._p1 *= radius
        self._p2 *= radius
        # Add current position
        self._p1 += self._c
        self._p2 += self._c
        # Set
        self.widget.SetPoint1(self._p1)
        self.widget.SetPoint2(self._p2)
        # Force the plane to be updated
#        self.widget.InvokeEvent("InteractionEvent")
        return
    
    def GetRadius(self):
        # Get into numpy vectors
        self.widget.GetCenter(self._c)
        self.widget.GetPoint1(self._p1)
        self.widget.GetPoint2(self._p2)

        # p1/2 now are relative to centre
        self._p1 -= self._c
        self._p2 -= self._c
        # Get norms
        r1 = N.sqrt(N.dot(self._p1, self._p1))
        r2 = N.sqrt(N.dot(self._p2, self._p2))
        return min(r1, r2)
    Radius = property(GetRadius, SetRadius)
#    @property
#    def IsRadiusValid(self):
#        radius = self.Radius
#        return N.isfinite(radius) and radius > 1e-6
    
#    @property
#    def CanShow(self):
#        return self.IsCentreValid and self.IsNormalValid and self.IsRadiusValid
    
    pass

class PlacedInlet(PlacedIolet):    
    def __init__(self):
        self.colour = (0.,1.,0.)
        PlacedIolet.__init__(self)
        return
    pass

class PlacedOutlet(PlacedIolet):    
    def __init__(self):
        self.colour = (1.,0.,0.)
        PlacedIolet.__init__(self)
        return
    pass

