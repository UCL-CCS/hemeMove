# This file is part of HemeLB and is Copyright (C)
# the HemeLB team and/or their institutions, as detailed in the
# file AUTHORS. This software is provided under the terms of the
# license in the file LICENSE.

from __future__ import absolute_import
import math
import wx

from .Layout import H, V, StretchSpacer, RectSpacer
from .VectorCtrl import VectorCtrl, VectorCtrlMapper
from .IoletListCtrl import IoletListCtrl

from ..Bindings.WxMappers import (
    WxWidgetMapper,
    WxWidgetEnabledMapper,
    NonObservingWxWidgetMapper,
    WxListCtrlMapper,
    WxListCtrlSelectionMapper,
)
from ..Bindings.Translators import (
    NoneToValueTranslator,
    FloatTranslator,
    QuickTranslator,
    Translator,
    Constraint,
)
from ..Bindings.Bindings import WxActionBinding

from ..Model.Profile import Profile


class ToolPanel(wx.Panel):
    """Tools Panel for the LHS of the window."""

    def __init__(self, controller, *args, **kwargs):
        wx.Panel.__init__(self, *args, **kwargs)
        self.controller = controller

        self.controlPanel = ControlPanel(controller, self)
        self.inputPanel = InputPanel(controller, self)
        self.ioletsPanel = IoletsPanel(controller, self)
        controller.BindValue(
            "HaveValidStlFile", WxWidgetEnabledMapper(self.ioletsPanel)
        )
        self.voxelPanel = VoxelPanel(controller, self)
        controller.BindValue("HaveValidStlFile", WxWidgetEnabledMapper(self.voxelPanel))
        self.seedPanel = SeedPanel(controller, self)
        controller.BindValue("HaveValidStlFile", WxWidgetEnabledMapper(self.seedPanel))
        self.outputPanel = OutputPanel(controller, self)
        controller.BindValue(
            "HaveValidStlFile", WxWidgetEnabledMapper(self.outputPanel)
        )

        layout = V(
            (self.controlPanel, 0, wx.EXPAND),
            (self.inputPanel, 0, wx.EXPAND),
            (self.ioletsPanel, 1, wx.EXPAND),
            (self.voxelPanel, 0, wx.EXPAND),
            (self.seedPanel, 0, wx.EXPAND),
            (self.outputPanel, 0, wx.EXPAND),
        )
        sizer = layout.create()
        self.SetSizer(sizer)
        return

    pass


class ControlPanel(wx.Panel):
    def __init__(self, controller, *args, **kwargs):
        """Setup the control panel."""
        wx.Panel.__init__(self, *args, **kwargs)
        self.controller = controller

        self.debugButton = wx.Button(self, label="DEBUG")
        controller.BindAction("Debug", WxActionBinding(self.debugButton, wx.EVT_BUTTON))

        self.openProfileButton = wx.Button(self, label="Open Profile")
        controller.BindAction(
            "LoadFromFile", WxActionBinding(self.openProfileButton, wx.EVT_BUTTON)
        )

        self.saveProfileButton = wx.Button(self, label="Save Profile")
        controller.BindAction(
            "ChooseSaveFile", WxActionBinding(self.saveProfileButton, wx.EVT_BUTTON)
        )

        self.generateButton = wx.Button(self, label="Generate")
        controller.BindValue(
            "IsReadyToGenerate", WxWidgetEnabledMapper(self.generateButton)
        )
        controller.BindAction(
            "Generate", WxActionBinding(self.generateButton, wx.EVT_BUTTON)
        )

        self.progressGauge = wx.Gauge(self)

        layout = V(
            (
                H(
                    self.debugButton,
                    StretchSpacer(),
                    self.openProfileButton,
                    self.saveProfileButton,
                    self.generateButton,
                    StretchSpacer(),
                ),
                0,
                wx.EXPAND,
            ),
            (self.progressGauge, 0, wx.EXPAND),
        )
        self.SetSizer(layout.create())

        return

    pass


class InputPanel(wx.Panel):
    def __init__(self, controller, *args, **kwargs):
        """Set up the input file part of the tool panel."""
        wx.Panel.__init__(self, *args, **kwargs)
        self.controller = controller

        iLabel = wx.StaticText(self, label="Input STL File")
        self.inputStlField = wx.TextCtrl(self, style=wx.TE_PROCESS_ENTER)

        uLabel = wx.StaticText(self, label="STL File Unit")
        self.inputStlUnitsChoice = wx.Choice(
            self, choices=[unit.Name for unit in Profile._UnitChoices]
        )

        controller.BindValue(
            "StlFile",
            WxWidgetMapper(
                self.inputStlField,
                "Value",
                wx.EVT_TEXT,
                translator=NoneToValueTranslator(""),
            ),
        )
        controller.BindValue(
            "HaveValidStlFile",
            NonObservingWxWidgetMapper(
                self.inputStlField,
                "BackgroundColour",
                translator=controller.validColourer,
            ),
        )

        self.inputStlButton = wx.Button(self, label="Choose")
        controller.BindAction(
            "ChooseStl", WxActionBinding(self.inputStlButton, wx.EVT_BUTTON)
        )
        controller.BindValue(
            "StlFileUnitId",
            WxWidgetMapper(self.inputStlUnitsChoice, "Selection", wx.EVT_CHOICE),
        )
        layout = V(
            iLabel,
            (
                H(
                    (V((self.inputStlField, 0, wx.EXPAND)), 1, wx.EXPAND),
                    self.inputStlButton,
                ),
                0,
                wx.EXPAND,
            ),
            H(uLabel, self.inputStlUnitsChoice),
        )
        self.SetSizer(layout.create())
        return

    pass


selectionToTrueTranslator = QuickTranslator(lambda x: x is not None, lambda x: None)


class LabelUnitTranslator(Translator):
    def __init__(self, template, inner=None):
        Translator.__init__(self, inner)
        self.template = template
        return

    def TranslateStage(self, val):
        return self.template % val.Abbrv

    def UntranslateState(self, val):
        return self.template % val.Abbrv

    pass


def MakeAndBindLengthUnitLabel(panel, controller, name):
    lbl = wx.StaticText(panel)
    controller.BindValue(
        "StlFileUnit",
        NonObservingWxWidgetMapper(
            lbl, "Label", translator=LabelUnitTranslator(name + " / %s")
        ),
    )
    return lbl


class IoletsDetailPanel(wx.Panel):
    def __init__(self, controller, *args, **kwargs):
        wx.Panel.__init__(self, *args, **kwargs)
        self.controller = controller

        nameLabel = wx.StaticText(self, label="Name")
        self.nameField = wx.TextCtrl(self)

        controller.BindValue(
            "Iolets.Selection.Name",
            WxWidgetMapper(
                self.nameField,
                "Value",
                wx.EVT_TEXT,
                translator=NoneToValueTranslator(""),
            ),
        )

        self.centreLabel = MakeAndBindLengthUnitLabel(self, controller, "Centre")

        self.centreVector = VectorCtrl(self)
        controller.BindValue(
            "Iolets.Selection.Centre",
            VectorCtrlMapper(self.centreVector, "Value", wx.EVT_TEXT),
        )

        self.radiusLabel = MakeAndBindLengthUnitLabel(self, controller, "Radius")
        self.radiusField = wx.TextCtrl(self)

        rad_trans = NoneToValueTranslator(
            math.nan,
            inner=Constraint(
                lambda x: math.isnan(x) or x > 0.0, inner=FloatTranslator()
            ),
        )
        controller.BindValue(
            "Iolets.Selection.Radius",
            WxWidgetMapper(
                self.radiusField, "Value", wx.EVT_TEXT, translator=rad_trans
            ),
        )

        self.placeButton = wx.Button(self, label="Place")
        controller.BindValue(
            "Pipeline.IoletPlaceButtonEnabled", WxWidgetEnabledMapper(self.placeButton)
        )
        controller.BindValue(
            "Pipeline.IoletPlaceButtonLabel",
            NonObservingWxWidgetMapper(self.placeButton, "Label"),
        )

        controller.BindAction(
            "Pipeline.IoletPlaceClicked",
            WxActionBinding(self.placeButton, wx.EVT_BUTTON),
        )

        normalLabel = wx.StaticText(self, label="Inward Normal")
        self.normalVector = VectorCtrl(self)
        controller.BindValue(
            "Iolets.Selection.Normal",
            VectorCtrlMapper(self.normalVector, "Value", wx.EVT_TEXT),
        )

        pressureLabel = wx.StaticText(self, label="Pressure / Pa")
        self.pressureVector = VectorCtrl(self)
        controller.BindValue(
            "Iolets.Selection.Pressure",
            VectorCtrlMapper(self.pressureVector, "Value", wx.EVT_TEXT),
        )

        self.pressureExpressionLabel = wx.StaticText(self, label="p =")
        controller.BindValue(
            "Iolets.Selection.PressureEquation",
            WxWidgetMapper(
                self.pressureExpressionLabel,
                "Label",
                wx.EVT_TEXT,
                translator=NoneToValueTranslator("p = "),
            ),
        )
        layout = V(
            (V(nameLabel, (self.nameField, 1, wx.EXPAND)), 0, wx.EXPAND),
            (self.placeButton, 0, wx.CENTRE),
            (V(self.radiusLabel, (self.radiusField, 1, wx.EXPAND)), 0, wx.EXPAND),
            (V(self.centreLabel, (self.centreVector, 1, wx.EXPAND)), 0, wx.EXPAND),
            (V(normalLabel, (self.normalVector, 1, wx.EXPAND)), 0, wx.EXPAND),
            (
                V(
                    pressureLabel,
                    (self.pressureVector, 1, wx.EXPAND),
                    self.pressureExpressionLabel,
                ),
                0,
                wx.EXPAND,
            ),
            RectSpacer(0, 2),
        )
        self.SetSizer(layout.create())

        controller.BindValue(
            "Iolets.SelectedIndex",
            WxWidgetEnabledMapper(self, translator=selectionToTrueTranslator),
        )

        return

    pass


class IoletsPanel(wx.Panel):
    def __init__(self, controller, *args, **kwargs):
        """Set up the inlets and outlets panel."""
        wx.Panel.__init__(self, *args, **kwargs)
        self.controller = controller

        ioletsLabel = wx.StaticText(self, label="Inlets and Outlets")

        self.ioletsListCtrl = IoletListCtrl(
            controller.Iolets, self, style=wx.LC_REPORT | wx.SUNKEN_BORDER
        )

        controller.BindValue("Iolets", WxListCtrlMapper(self.ioletsListCtrl))

        controller.BindValue(
            "Iolets.SelectedIndex", WxListCtrlSelectionMapper(self.ioletsListCtrl)
        )

        self.addInletButton = wx.Button(self, label="Add Inlet")
        controller.BindAction(
            "Iolets.AddInlet", WxActionBinding(self.addInletButton, wx.EVT_BUTTON)
        )

        self.addOutletButton = wx.Button(self, label="Add Outlet")
        controller.BindAction(
            "Iolets.AddOutlet", WxActionBinding(self.addOutletButton, wx.EVT_BUTTON)
        )

        self.removeIoletButton = wx.Button(self, label="Remove")
        controller.BindValue(
            "Iolets.SelectedIndex",
            WxWidgetEnabledMapper(
                self.removeIoletButton, translator=selectionToTrueTranslator
            ),
        )
        controller.BindAction(
            "Iolets.RemoveIolet", WxActionBinding(self.removeIoletButton, wx.EVT_BUTTON)
        )

        self.detail = IoletsDetailPanel(controller, self)

        layout = V(
            ioletsLabel,
            (
                H(
                    (V((self.ioletsListCtrl, 1, wx.EXPAND)), 1, wx.EXPAND),
                    (V(StretchSpacer(), (self.detail, 0, wx.EXPAND)), 2, wx.EXPAND),
                ),
                1,
                wx.EXPAND,
            ),
            (
                H(self.addInletButton, self.addOutletButton, self.removeIoletButton),
                0,
                wx.EXPAND,
            ),
            RectSpacer(0, 1),
        )
        self.SetSizer(layout.create())
        return

    pass


class VoxelPanel(wx.Panel):
    def __init__(self, controller, *args, **kwargs):
        """Setup up the voxel size panel."""
        wx.Panel.__init__(self, *args, **kwargs)
        self.controller = controller

        self.voxelLabel = MakeAndBindLengthUnitLabel(self, controller, "Voxel size")
        #        controller.BindValue('StlFileUnitId',
        #                              WxWidgetMapper(voxelLabel,
        #                                             'Label', wx.EVT_TEXT,
        #                                             translator=QuickTranslator(lambda i: 'Voxel size (%s)' % Profile._UnitChoices[i].Abbrv,
        #                                                                        lambda x: 1)))
        self.voxelSizeField = wx.TextCtrl(self)
        self.controller.BindValue(
            "VoxelSize",
            WxWidgetMapper(
                self.voxelSizeField, "Value", wx.EVT_TEXT, translator=FloatTranslator()
            ),
        )

        self.voxelResetButton = wx.Button(self, label="Reset")
        self.controller.BindAction(
            "ResetVoxelSize", WxActionBinding(self.voxelResetButton, wx.EVT_BUTTON)
        )

        layout = V(
            self.voxelLabel,
            (
                H(
                    (V((self.voxelSizeField, 0, wx.EXPAND)), 1, wx.EXPAND),
                    self.voxelResetButton,
                ),
                1,
                wx.EXPAND,
            ),
        )
        #        layout = V(
        #            xmlLabel,
        #            (H((V((self.xmlField, 0, wx.EXPAND)), 1, wx.EXPAND),
        #                 self.xmlChooseButton
        #                 ), 0, wx.EXPAND)
        #            )
        self.SetSizer(layout.create())
        return

    pass


class SeedPanel(wx.Panel):
    def __init__(self, controller, *args, **kwargs):
        """Setup the seed point panel."""
        wx.Panel.__init__(self, *args, **kwargs)
        self.controller = controller

        self.seedLabel = MakeAndBindLengthUnitLabel(self, controller, "Seed position")
        self.seedVector = VectorCtrl(self)

        controller.BindValue(
            "SeedPoint", VectorCtrlMapper(self.seedVector, "Value", wx.EVT_TEXT)
        )
        controller.BindValue(
            "HaveValidSeedPoint",
            NonObservingWxWidgetMapper(
                self.seedVector, "BackgroundColour", translator=controller.validColourer
            ),
        )

        self.seedPlaceButton = wx.Button(self, label="Place")
        controller.BindValue(
            "Pipeline.SeedPlaceButtonEnabled",
            WxWidgetEnabledMapper(self.seedPlaceButton),
        )
        controller.BindValue(
            "Pipeline.SeedPlaceButtonLabel",
            NonObservingWxWidgetMapper(self.seedPlaceButton, "Label"),
        )

        controller.BindAction(
            "Pipeline.SeedPlaceClicked",
            WxActionBinding(self.seedPlaceButton, wx.EVT_BUTTON),
        )

        layout = V(
            self.seedLabel,
            (H((self.seedVector, 1, wx.EXPAND), self.seedPlaceButton), 1, wx.EXPAND),
        )
        self.SetSizer(layout.create())
        return

    pass


class OutputPanel(wx.Panel):
    def __init__(self, controller, *args, **kwargs):
        """Setup the output file panel."""
        wx.Panel.__init__(self, *args, **kwargs)
        self.controller = controller

        geometryLabel = wx.StaticText(self, label="Output geometry")
        self.geometryField = wx.TextCtrl(self)
        controller.BindValue(
            "OutputGeometryFile",
            WxWidgetMapper(
                self.geometryField,
                "Value",
                wx.EVT_TEXT,
                translator=NoneToValueTranslator(""),
            ),
        )
        controller.BindValue(
            "HaveValidOutputGeometryFile",
            NonObservingWxWidgetMapper(
                self.geometryField,
                "BackgroundColour",
                translator=controller.validColourer,
            ),
        )

        self.geometryChooseButton = wx.Button(self, label="Choose")
        controller.BindAction(
            "ChooseOutputGeometryFile",
            WxActionBinding(self.geometryChooseButton, wx.EVT_BUTTON),
        )

        xmlLabel = wx.StaticText(self, label="Output xml")
        self.xmlField = wx.TextCtrl(self)
        controller.BindValue(
            "OutputXmlFile",
            WxWidgetMapper(
                self.xmlField,
                "Value",
                wx.EVT_TEXT,
                translator=NoneToValueTranslator(""),
            ),
        )
        controller.BindValue(
            "HaveValidOutputXmlFile",
            NonObservingWxWidgetMapper(
                self.xmlField, "BackgroundColour", translator=controller.validColourer
            ),
        )

        self.xmlChooseButton = wx.Button(self, label="Choose")
        controller.BindAction(
            "ChooseOutputXmlFile", WxActionBinding(self.xmlChooseButton, wx.EVT_BUTTON)
        )

        layout = V(
            xmlLabel,
            (
                H(
                    (V((self.xmlField, 0, wx.EXPAND)), 1, wx.EXPAND),
                    self.xmlChooseButton,
                ),
                0,
                wx.EXPAND,
            ),
            geometryLabel,
            (
                H(
                    (V((self.geometryField, 0, wx.EXPAND)), 1, wx.EXPAND),
                    self.geometryChooseButton,
                ),
                0,
                wx.EXPAND,
            ),
        )

        self.SetSizer(layout.create())
        return

    pass
