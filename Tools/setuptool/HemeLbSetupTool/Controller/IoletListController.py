import wx

from HemeLbSetupTool.Bindings.ListController import ListController, HasListKeys
from HemeLbSetupTool.Controller.VectorController import VectorController
from HemeLbSetupTool.Controller.IoletController import IoletController
from HemeLbSetupTool.Model.Iolets import Inlet, Outlet
import pdb
class IoletListController(ListController):
    def __init__(self, delegate):
        ListController.__init__(self, delegate, SelectionControllerClass=IoletController)
        self.nInlets = 0
        self.nOutlets = 0
        return
    
    def AddInlet(self):
        self.nInlets += 1
        self.delegate.append(Inlet(Name='Inlet%d' % (self.nInlets)))
        self.SelectedIndex = len(self.delegate)-1
        return
    
    def AddOutlet(self):
        self.nOutlets += 1
        self.delegate.append(Outlet(Name='Outlet%d' % (self.nOutlets)))
        self.SelectedIndex = len(self.delegate)-1
        return

    def RemoveIolet(self):
        pdb.set_trace()
        if self.SelectedIndex is None:
            return
        del self.delegate[self.SelectedIndex]
        return
    
    pass

class HasIoletListKeys(HasListKeys):
    """Mixin for ObjectController subclasses with IoletList keys.
    """
    BindMethodDispatchTable = ((IoletListController, 'BindList'),)

    # def BindList(self, modelKey, widgetMapper):
    #     """We need to bind the selection and deal with add/remove/update.
    #     """
    #     self.BindComplexValue(modelKey, ListContentsMapper, (),
    #                           ValueBinding, widgetMapper)
        
    #     return

    def DefineIoletListKey(self, name):
        """Typically used in the subclass __init__ method to easily
        mark a key as being a List and hence needing a ListController
        to manage it.
        """
        setattr(self, name,
                IoletListController(getattr(self.delegate, name))
                )
        return
    
    pass
