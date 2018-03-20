#! /usr/bin/env python
# -*- coding: utf-8 -*-

"""
The CroCo Cross-Link Converter GUI

Graphical interface to convert results from data analysis of chemical cross-linking /
mass-spectrometry experiments.

This script creates the GUI in wxPython (https://wxpython.org/pages/overview/)
"""
import sys, os

croco_script = sys.argv[0]
# check if programme was called via symlink
if os.path.islink(croco_script):
    croco_script = os.readlink(croco_script)
# dir is the directory above the bin-dir
croco_dir = os.path.abspath(os.path.join(os.path.dirname(croco_script), '..'))

sys.path.append(os.path.abspath(os.path.join('..', croco_dir)))

import wx

import croco
from pandas import read_csv

class CrocoMainWindow(wx.Frame):

    def __init__(self):
        # ensure the parent's __init__ is called
        wx.Frame.__init__(self,
                          None,
                          wx.ID_ANY,
                          title='The CroCo cross-link converter')

        # create a panel in the frame
        self.panel = wx.Panel(self)

        # define the croco read and write options and link to the modules
        self.availReads = {'pLink1': 'croco.pLink1.Read',
                           'pLink2': 'croco.pLink2.Read',
                           'Kojak': 'croco.Kojak.Read',
                           'xQuest': 'croco.xQuest.Read',
                           'xTable': 'read_csv'}
        self.availWrites =  {'xTable': 'croco.xTable.Write',
                             'xVis': 'croco.xVis.Write',
                             'xiNet': 'croco.xiNET.Write',
                             'DynamXL': 'croco.DynamXL.Write',
                             'xWalk': 'croco.xWalk.Write'}

        # set triggers allowing the start button to unhide
        self.readSet = False
        self.writeSet = False

        # create a menu bar
        self.makeMenuBar()

        # and a status bar
        self.CreateStatusBar()
        self.SetStatusText("Welcome to CroCo!")

        # load widgets
        # ALWAYS AFTER LOADING MENU AND STATUS BAR!!!!
        self.createWidgets()

    def createWidgets(self):

        ## define the widgets
        # define the contents of the first sizer
        self.inputButton = wx.Button(self.panel, label='Load file(s)')
        self.inputButton.Enable(False)
        self.readFormat = wx.Choice(self.panel, choices=list(self.availReads.keys()))
        self.outputButton = wx.Button(self.panel, label='Write to')
        self.outputButton.Enable(False)
        self.writeFormat = wx.Choice(self.panel, choices=list(self.availWrites.keys()))
        
        self.controlStart = wx.Button(self.panel, label='Start')
        self.controlStart.Enable(False)
        controlQuit = wx.Button(self.panel, label='Quit')

        gauge = wx.Gauge(self.panel, range=100)

        ## create Connects

        self.readFormat.Bind(wx.EVT_CHOICE, self.OnReadFormat)
        self.writeFormat.Bind(wx.EVT_CHOICE, self.OnWriteFormat)

        self.inputButton.Bind(wx.EVT_BUTTON, self.OnOpenSwitch)
        self.outputButton.Bind(wx.EVT_BUTTON, self.OnOutputDir)

        controlQuit.Bind(wx.EVT_BUTTON, self.OnExit)

        ## define the layout
        # define the sizers (main widget layouts) used in the app
        topSizer = wx.BoxSizer(wx.VERTICAL)
        loadSizer = wx.FlexGridSizer(rows=2, cols=2, vgap=10, hgap=10)
        loadSizerContainer = wx.BoxSizer(wx.HORIZONTAL)
        controlSizer = wx.BoxSizer(wx.HORIZONTAL)

        # assign the widgets to the sizers
        loadSizer.AddMany([(self.readFormat, 1, wx.EXPAND),
                           self.inputButton,
                           (self.writeFormat, 1, wx.EXPAND),
                           self.outputButton])

        # allow to resize the second row
        loadSizer.AddGrowableCol(1, 1)
        # contain the grid in a BoxSizer to set borders individually
        loadSizerContainer.Add(loadSizer, proportion=1, flag=wx.ALL|wx.EXPAND)

        # Set start and quit buttons
        controlSizer.Add(self.controlStart, 1, wx.RIGHT, 5)
        controlSizer.Add(controlQuit, 0, wx.LEFT | wx.ALIGN_RIGHT, 5)

        # assign lower sizers to the top-level sizer
        topSizer.Add(loadSizerContainer, 0, wx.ALL|wx.EXPAND, 5)
        # topSizer.Add(reviewSizer, 0, wx.ALL|wx.EXPAND, 5)
        topSizer.Add(wx.StaticLine(self.panel), 0, wx.ALL|wx.EXPAND, 5)
        topSizer.Add(controlSizer, 0, wx.ALL|wx.EXPAND, 5)
        topSizer.Add(gauge, 0,wx.ALL|wx.EXPAND, 5)

        # assign the top sizer to the panel i.e. main layout instance
        self.panel.SetSizer(topSizer)
        topSizer.Fit(self)
        
    def makeMenuBar(self):
        """
        A menu bar is composed of menus, which are composed of menu items.
        This method builds a set of menus and binds handlers to be called
        when the menu item is selected.
        """

        # Make a file menu
        fileMenu = wx.Menu()
        
        # The "\t..." syntax defines an accelerator key that also triggers
        # the same event
        menu_load = fileMenu.Append(-1, "&Load file\tCtrl-O",
                "Load a new file")
        fileMenu.AppendSeparator()

        # When using a stock ID we don't need to specify the menu item's
        # label
        exitItem = fileMenu.Append(wx.ID_EXIT)

        # Now a help menu for the about item
        helpMenu = wx.Menu()
        aboutItem = helpMenu.Append(wx.ID_ABOUT)

        # Make the menu bar and add the two menus to it. The '&' defines
        # that the next letter is the "mnemonic" for the menu item. On the
        # platforms that support it those letters are underlined and can be
        # triggered from the keyboard.
        menuBar = wx.MenuBar()
        menuBar.Append(fileMenu, "File")
        menuBar.Append(helpMenu, "Help")

        # Give the menu bar to the frame
        self.SetMenuBar(menuBar)

        # Finally, associate a handler function with the EVT_MENU event for
        # each of the menu items. That means that when that menu item is
        # activated then the associated handler function will be called.
        self.Bind(wx.EVT_MENU, self.OnLoad, menu_load)
        self.Bind(wx.EVT_MENU, self.OnExit,  exitItem)
        self.Bind(wx.EVT_MENU, self.OnAbout, aboutItem)

    ## Functions

    def OnReadFormat(self, event):
        self.theReadFormat = self.readFormat.GetString(self.readFormat.GetSelection())
        print('Reading {} format'.\
            format(self.theReadFormat))
        self.inputButton.Enable(True)
        
    def OnWriteFormat(self, event):
        self.theWriteFormat = self.writeFormat.GetString(self.writeFormat.GetSelection())
        print('Writing {} format'.\
            format(self.theWriteFormat))
        self.outputButton.Enable(True)

    def OnOpenSwitch(self, event):
        if 'pLink' in self.theReadFormat:
            self.OnOpenDir()
        else:
            self.OnOpenFile()

    def OnOpenFile(self):
       dlg = wx.FileDialog(self, 
                           message="Choose a file",
                           defaultDir=os.getcwd(),
                           defaultFile="",
                           wildcard="*.*",
                           style=wx.FD_OPEN)
       if dlg.ShowModal() == wx.ID_OK:
            path = dlg.GetPath()
            self.inpath = os.path.basename(path)
       dlg.Destroy()
       
       self.readSet = True
       if self.writeSet:
           self.controlStart.Enable(True)

    def OnOpenDir(self):
        dlg = wx.DirDialog(self, "Choose a directory:", style=wx.DD_DEFAULT_STYLE | wx.DD_NEW_DIR_BUTTON)
        if dlg.ShowModal() == wx.ID_OK:
            path = dlg.GetPath()
            self.inpath = os.path.basename(path)
        dlg.Destroy()
        
        self.readSet = True
        if self.writeSet:
            self.controlStart.Enable(True)
       
    def OnOutputDir(self, event):
        dlg = wx.DirDialog(self, "Choose a directory:", style=wx.DD_DEFAULT_STYLE | wx.DD_NEW_DIR_BUTTON)
        if dlg.ShowModal() == wx.ID_OK:
            path = dlg.GetPath()
            self.outpath = os.path.basename(path)
        dlg.Destroy()
        
        self.writeSet = True
        if self.readSet:
            self.controlStart.Enable(True)

    def OnExit(self, event):
        """Close the frame, terminating the application."""
        dlg = wx.MessageDialog(self, 
            "Do you really want to close this application?",
            "Confirm Exit", wx.OK|wx.CANCEL|wx.ICON_QUESTION)
        result = dlg.ShowModal()
        dlg.Destroy()
        if result == wx.ID_OK:
            self.Destroy()

    def OnLoad(self, event):
        pass

    def onCancel(self, event):
        self.closeProgram()

    def OnAbout(self, event):
        """Display an About Dialog"""
        wx.MessageBox("This is a wxPython Hello World sample",
                      "About Hello World 2",
                      wx.OK|wx.ICON_INFORMATION)


if __name__ == '__main__':
    # When this module is run (not imported) then create the app, the
    # frame, show it, and start the event loop.
    app = wx.App()
    frm = CrocoMainWindow().Show()
    app.MainLoop()