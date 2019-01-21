#! /usr/bin/env python
# -*- coding: utf-8 -*-

"""
The CroCo Cross-Link Converter GUI

Graphical interface to convert results from data analysis of chemical cross-linking /
mass-spectrometry experiments.

This script creates the GUI in wxPython (https://wxpython.org/pages/overview/)
"""
import os, sys

import wx
import wx.adv

import croco
#from pandas import read_csv

class CroCoMainFrame(wx.Frame):

    def __init__(self):
        # ensure the parent's __init__ is called
        wx.Frame.__init__(self,
                          None,
                          wx.ID_ANY,
                          style= wx.SYSTEM_MENU | wx.CAPTION | wx.CLOSE_BOX,
                          title='The CroCo cross-link converter')

        # create a panel in the frame
        self.panel = wx.Panel(self)

        ## setting the icon for the frame
        # in case of calling croco from the source folder structure...
        file_dir, file_name = os.path.split(__file__)
        if os.path.exists(os.path.join(file_dir,
                                       './croco/data/croco_logo.ico')):
            iconFile = os.path.abspath(os.path.join(file_dir,
                                                      './croco/data/croco_logo.ico'))
        # ... or calling from a exe-file in a folder-setup with the data folder at top-level
        elif os.path.exists('./croco/data/croco_logo.ico'):
            iconFile = os.path.abspath('./croco/data/croco_logo.ico')
        # ... or calling from within a single bundled exe-file
        else:
            try:
                # PyInstaller creates a temp folder and stores its path in _MEIPASS
                base_path = sys._MEIPASS
                iconFile =  os.path.abspath(\
                    os.path.join(base_path, './data/croco_logo.ico'))
            # ... or something went wrong
            except:
                raise Exception('croco_logo.ico not found')
        ## end setting the icon

        icon = wx.Icon(iconFile, wx.BITMAP_TYPE_ICO)
        self.SetIcon(icon)

        # set starting path for folder selection GUI
        self.currentPath = os.getcwd()

        # define the croco read and write options and link to the modules
        # each value of a dict is a list of [options, functionToCall]
        self.availReads = {'pLink1': [croco.pLink1.Read, []],
                           'pLink2': [croco.pLink2.Read, []],
                           'Kojak': [croco.Kojak.Read, [('Kojak rawfile',
                                                         'input',
                                                         'Please provide the name of the rawfile for the data')]],
                           'Kojak + Percolator': [croco.KojakPercolator.Read, []],
                           'StavroX': [croco.StavroX.Read, [('SSF file',
                                                             'file',
                                                             'Please provide an SSF file as returned by StavroX')]],
                           'Xi': [croco.Xi.Read, []],
                           'Xi + XiFDR': [croco.XiSearchFDR.Read,[('Xi search file',
                                                                   'file',
                                                                   'Please provide Xi output file')]],
                           'xQuest': [croco.xQuest.Read, []],
                           'xTable': [croco.xTable.Read, []]}

        """
        availWrites is a dict with lists as values for every possible writing
        operation as key. The lists are structured as follows
        1: function to call
        2: options for displaying a related widget (is a list of tuples)
            1: Label for the widget
            2. type of the widget (file/dir/input/check)
            3. Help text
            4. Optional multu-purpose field: e.g. default value for input fields
        """
        self.availWrites =  {'xTable': [croco.xTable.Write, []],
                             'xVis': [croco.xVis.Write, []],
                             'xiNet': [croco.xiNET.Write, []],
                             'DynamXL': [croco.DynamXL.Write, []],
                             'xWalk': [croco.xWalk.Write, [('PDB to map xlinks to',
                                                            'file',
                                                            'Please provide a PDB file you want to map the xlinks to.\n'+\
                                                            'No special format is required as only the filename is used'),
                                                           ('Offset',
                                                            'input',
                                                            'Difference between the index in the PDB file and the xTable ' +\
                                                            '(e.g. 1 for 12 in the xTable and 13 in the PDB)'),
                                                           ('Chains',
                                                            'input',
                                                            'Map chains to protein names (e.g. Protein1:AB, Protein2:C). '+\
                                                            'Use exactly the same naming scheme as the xTable file'),
                                                           ('PDB Atom code',
                                                            'input',
                                                            'Provide an PDB atom code (e.g. CB) for distance calculation')]],
                             'pLabel': [croco.pLabel.Write, [('Dir containing mgf files',
                                                              'dir',
                                                              'Provide path to a folder containing mgf-files corresponding to the'+\
                                                              'rawfile names in the input file'),
                                                             ('Xlinker as referenced by pLabel',
                                                              'input',
                                                              'Provide a cross-linker name (e.g. BS3) that is used by pLabel'+\
                                                              'to calculate the potentially linked amino acids',
                                                              'BS3'),
                                                             ('Merge mgf files',
                                                              'check',
                                                              'Tick to merge mgf files into one file containign' +
                                                              'only the spectra  mentioned in the (merged) pLabel')]],
                              'customTable': [croco.customTable.Write, [('Custom Template',
                                                                         'file',
                                                                         'Provide tample file for parsing')]]}

        # defines the column headers required for xtable output
        self.col_order = ['rawfile', 'scanno', 'prec_ch',
                          'pepseq1', 'xlink1',
                          'pepseq2', 'xlink2', 'xtype',
                          'modmass1', 'modpos1', 'mod1',
                          'modmass2', 'modpos2', 'mod2',
                          'prot1', 'xpos1', 'prot2',
                          'xpos2', 'type', 'score', 'ID', 'pos1', 'pos2', 'decoy']

        # set triggers allowing the start button to unhide
        self.readSet = False
        self.writeSet = False

        # create a menu bar
        self.makeMenuBar()

        # load widgets
        # ALWAYS AFTER LOADING MENU AND STATUS BAR!!!!
        self.createWidgets()

        # arguments to pass at function call
        self.inputArgs = ''
        self.outputArgs = ''

    def createWidgets(self):

        ## define the widgets

        # define the contents of the first sizer

        input_lbl = wx.StaticText(self.panel,wx.ID_ANY, label='Input', style=wx.ALIGN_CENTER )
        self.inputButton = wx.Button(self.panel, label='Load file(s)')
        self.inputButton.Enable(False)
        self.readFormat = wx.Choice(self.panel, choices=sorted(list(self.availReads.keys())))

        output_lbl = wx.StaticText(self.panel,wx.ID_ANY, label='Output', style=wx.ALIGN_CENTER )
        self.outputButton = wx.Button(self.panel, label='Write to')
        self.outputButton.Enable(False)
        self.writeFormat = wx.Choice(self.panel, choices=sorted(list(self.availWrites.keys())))

        self.compactTableCheck = wx.CheckBox(self.panel, label='Compact xTable')

        self.controlStart = wx.Button(self.panel, label='Start')
        self.controlStart.Enable(False)
        controlQuit = wx.Button(self.panel, label='Quit')

        helpButton = wx.ContextHelpButton(self.panel)

        # progress bar
        # gauge = wx.Gauge(self.panel, range=100)

        ## create Connects

        self.readFormat.Bind(wx.EVT_CHOICE, self.OnReadFormat)
        self.readFormat.Bind(wx.EVT_HELP,
                             lambda evt: self.Info('Select data-format for input', caption='Help'))
        self.writeFormat.Bind(wx.EVT_CHOICE, self.OnWriteFormat)
        self.writeFormat.Bind(wx.EVT_HELP,
                              lambda evt: self.Info('Select data-format for output', caption='Help'))

        self.inputButton.Bind(wx.EVT_BUTTON, self.OnOpenSwitch)
        self.inputButton.Bind(wx.EVT_HELP,
                              lambda evt: self.Info('Opens a dialog to select file for input', caption='Help'))
        self.outputButton.Bind(wx.EVT_BUTTON, self.onOutputDir)
        self.outputButton.Bind(wx.EVT_HELP,
                              lambda evt: self.Info('Opens a dialog to select output dir', caption='Help'))

        controlQuit.Bind(wx.EVT_BUTTON, self.onExit)
        self.controlStart.Bind(wx.EVT_BUTTON, self.onStart)
        self.controlStart.Bind(wx.EVT_HELP,
                               lambda evt: self.Info('Start the conversion', caption='Help'))
        ## define the layout
        # define the sizers (main widget layouts) used in the app
        topSizer = wx.BoxSizer(wx.VERTICAL)
        loadSizer = wx.FlexGridSizer(rows=2, cols=3, vgap=10, hgap=10)
        loadSizerContainer = wx.BoxSizer(wx.HORIZONTAL)
        controlSizer = wx.BoxSizer(wx.HORIZONTAL)
        checkSizer = wx.BoxSizer(wx.HORIZONTAL)

        # assign the widgets to the sizers
        loadSizer.AddMany([input_lbl,
                          (self.readFormat, 1, wx.EXPAND),
                          self.inputButton,
                          output_lbl,
                          (self.writeFormat, 1, wx.EXPAND),
                          self.outputButton])

        # allow to resize the second row
        loadSizer.AddGrowableCol(1, 1)
        # contain the grid in a BoxSizer to set borders individually
        loadSizerContainer.Add(loadSizer, proportion=1, flag=wx.ALL|wx.EXPAND)

        # Set checkboxes
        checkSizer.Add(self.compactTableCheck, 1, wx.ALL|wx.EXPAND, 5)

        # Set start and quit buttons
        controlSizer.Add(self.controlStart, 1, wx.RIGHT, 5)
        controlSizer.Add(controlQuit, 0, wx.LEFT | wx.ALIGN_RIGHT, 5)
        controlSizer.Add(helpButton, 0, wx.LEFT | wx.ALL | wx.EXPAND, 5)

        # assign lower sizers to the top-level sizer
        topSizer.Add(loadSizerContainer, 0, wx.ALL|wx.EXPAND, 5)
        # topSizer.Add(reviewSizer, 0, wx.ALL|wx.EXPAND, 5)
        topSizer.Add(checkSizer, 0, wx.EXPAND|wx.RIGHT, 5)
        topSizer.Add(wx.StaticLine(self.panel), 0, wx.ALL|wx.EXPAND, 5)
        topSizer.Add(controlSizer, 0, wx.EXPAND|wx.RIGHT, 5)
        # topSizer.Add(gauge, 0,wx.ALL|wx.EXPAND, 5)

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
        # menu_load = fileMenu.Append(-1, "&Load file\tCtrl-O",
        #         "Load a new file")
        # fileMenu.AppendSeparator()

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
        # self.Bind(wx.EVT_MENU, self.OnLoad, menu_load)
        self.Bind(wx.EVT_MENU, self.onExit,  exitItem)
        self.Bind(wx.EVT_MENU, self.onAbout, aboutItem)

    ## GUI Functions

    def OnReadFormat(self, event):
        self.theReadFormat = self.readFormat.GetString(self.readFormat.GetSelection())
        print('[onReadFormat] Reading {} format'.\
            format(self.theReadFormat))
        self.inputButton.Enable(True)

    def OnWriteFormat(self, event):
        self.theWriteFormat = self.writeFormat.GetString(self.writeFormat.GetSelection())
        print('[onWriteFormat] Writing {} format'.\
            format(self.theWriteFormat))
        self.outputButton.Enable(True)

    def OnOpenSwitch(self, event):
        if 'pLink' in self.theReadFormat:
            self.onInputDir()
        else:
            self.onInputFile()

    def onInputFile(self):
       dlg = wx.FileDialog(self,
                           message="Choose one or multiple files for input",
                           defaultDir=self.currentPath,
                           defaultFile="",
                           wildcard="*.*",
                           style=wx.FD_MULTIPLE)
       if dlg.ShowModal() == wx.ID_OK:
            self.theInput = dlg.GetPaths()
            self.currentPath = os.path.dirname(self.theInput[0])
            print('[onInputFile] Loaded {}'.format(', '.join(self.theInput)))
       dlg.Destroy()

       self.readSet = True
       if self.writeSet:
           self.controlStart.Enable(True)

    def onInputDir(self):
        dlg = wx.DirDialog(self,
                           message="Choose one or multiple directories for Input:",
                           defaultPath=self.currentPath,
                           style=wx.DD_DEFAULT_STYLE | wx.DD_DIR_MUST_EXIST)

        if dlg.ShowModal() == wx.ID_OK:
            # self.theInput is always a list of paths
            self.theInput = [dlg.GetPath()]
            self.currentPath = dlg.GetPath()
            print('[onInputDir] Loaded {}'.format(self.currentPath))

        dlg.Destroy()

        self.readSet = True
        if self.writeSet:
            self.controlStart.Enable(True)

    def onOutputDir(self, event):
        dlg = wx.DirDialog(self,
                           message="Choose a directory:",
                           defaultPath=self.currentPath,
                           style=wx.DD_DEFAULT_STYLE | wx.DD_NEW_DIR_BUTTON)
        if dlg.ShowModal() == wx.ID_OK:
            self.currentPath = self.theOutput = dlg.GetPath()
            print('[onOutputDir] Loaded {}'.format(self.currentPath))

        dlg.Destroy()

        self.writeSet = True
        if self.readSet:
            self.controlStart.Enable(True)

    def onExit(self, event):
        """Close the frame, terminating the application."""
        dlg = wx.MessageDialog(self,
            "Do you really want to close this application?",
            "Confirm Exit", wx.OK|wx.CANCEL|wx.ICON_QUESTION)
        result = dlg.ShowModal()
        dlg.Destroy()
        if result == wx.ID_OK:
            self.Destroy()

    def onCancel(self, event):
        self.closeProgram()

    def onAbout(self, event):

        aboutInfo = wx.adv.AboutDialogInfo()
        aboutInfo.SetName("The CroCo cross-link converter")
        aboutInfo.SetVersion('0.5.11')
        aboutInfo.SetDescription("Graphical interface to convert results from "+\
                                 "data analysis of chemical cross-linking "+\
                                 "mass-spectrometry experiments.")
        aboutInfo.SetCopyright("(C) 2018")
        aboutInfo.SetWebSite("www.halomem.de")
        aboutInfo.AddDeveloper("Julian Bender (jub@halomem.de)")

        wx.adv.AboutBox(aboutInfo)

    def Warning(self, message, caption = 'Warning!'):
        print('[WARN] {}'.format(message))
        dlg = wx.MessageDialog(self, message, caption, wx.OK | wx.ICON_WARNING)
        dlg.ShowModal()
        dlg.Destroy()

    def Info(self, message, caption = 'CroCo'):
        dlg = wx.MessageDialog(self, message, caption, wx.OK | wx.ICON_INFORMATION)
        dlg.ShowModal()
        dlg.Destroy()

    ## Procedural Functions

    def onStart(self, event):
        """
        Show dialog if additional user input is required.
        Otherwise start the conversion
        """
        # Check if there are options to ask for
        if len(self.availReads[self.theReadFormat][1]) > 0 or\
            len(self.availWrites[self.theWriteFormat][1]) > 0:

            # init the OptionsWindow as child
            OptionsFrame = CroCoOptionsFrame(self)

            # set the variables in the child window
            OptionsFrame.InputOptionsToAsk = self.availReads[self.theReadFormat][1]
            OptionsFrame.OutputOptionsToAsk = self.availWrites[self.theWriteFormat][1]

            # update the controls according to the variables
            OptionsFrame.updateOptions()
            # show the window
            OptionsFrame.Show()

        else:
            # reset args lists
            self.inputArgs = ''
            self.outputArgs = ''

            self.onRun(event)

    def onRun(self, event):
        """
        Collect all necessary information from self and start the
        conversion by calling the actual conversion script
        """

        print('[onRun] Going to convert {} from {} '.format(', '.join(self.theInput),
                                                   self.theReadFormat) +
               'format to {} format'.format(self.theWriteFormat))

        was_error = False

        # Displays a busy cursor during the run of the programme
        wait = wx.BusyCursor()

        for f in self.theInput:
            try:
                if len(self.inputArgs) > 0 :
                    print('[onRun] Found input args: "{}"'.format(self.inputArgs))
                    xtable = self.availReads[self.theReadFormat][0](f, *self.inputArgs)
                else:
                    print('[onRun] No extra input arguments provided.')
                    xtable = self.availReads[self.theReadFormat][0](f)
                print('[onRun] Table successfully read: {}'.format(f))
            except Exception as e:
                self.Warning('Error while reading Input-file: ' + str(e))

            print('[onRun] xTable read from input: {}'.format(', '.join(xtable.columns)))

            # Compact the xTable if checkbox is checked
            xtable = croco.HelperFunctions.applyColOrder(xtable,
                                                         col_order=self.col_order,
                                                         compact=self.compactTableCheck.GetValue())

            # if no user-defined output dir use current
            if self.theOutput == '':
                self.theOutput = os.path.dirname(f)

            # set filename for output file
            fname = os.path.splitext(os.path.split(f)[1])[0] + '_' + self.theReadFormat +\
                    '_to_' + self.theWriteFormat

            # generate output path w/o extension
            outpath = os.path.join(self.theOutput, fname)

            try:
                print('[onRun] Writing table in {} format to {}'.format(self.theWriteFormat, outpath))
                if len(self.outputArgs) > 0:
                    print('[onRun] Found output args: "{}"'.format(self.outputArgs))
                    self.availWrites[self.theWriteFormat][0](xtable, outpath, *self.outputArgs)
                else:
                    print('[onRun] No extra input arguments provided.')
                    self.availWrites[self.theWriteFormat][0](xtable, outpath)
                print('[onRun] Table successfully written!')
            except Exception as e:
                self.Warning('[onRun] Conversion of {} was '.format(f) +
                                   'not successfull:{}'.format(str(e)))
                was_error = True
                break
        # ends busy cursor
        del wait

        if not was_error:
            self.Info('File(s) successfully written ' +
                     'to {}!'.format(outpath),
                     caption='Success!')

class CroCoOptionsFrame(wx.Frame):
    """
    This is a Child Frame to CroCoMainWindow asking the user for input
    regarding the possible optinos of a submodule
    """

    def __init__(self, parent):
        wx.Frame.__init__(self,
                          None,
                          wx.ID_ANY,
                          style= wx.SYSTEM_MENU | wx.CAPTION | wx.CLOSE_BOX,
                          title='Options')

        self.panel = wx.Panel(self)

        ## setting the icon for the frame
        # in case of calling croco from the source folder structure...
        file_dir, file_name = os.path.split(__file__)
        if os.path.exists(os.path.join(file_dir,
                                       './croco/data/croco_logo.ico')):
            iconFile = os.path.abspath(os.path.join(file_dir,
                                                      './croco/data/croco_logo.ico'))
        # ... or calling from a exe-file in a folder-setup with the data folder at top-level
        elif os.path.exists('./croco/data/croco_logo.ico'):
            iconFile = os.path.abspath('./croco/data/croco_logo.ico')
        # ... or calling from within a single bundled exe-file
        else:
            try:
                # PyInstaller creates a temp folder and stores its path in _MEIPASS
                base_path = sys._MEIPASS
                iconFile =  os.path.abspath(\
                    os.path.join(base_path, './data/croco_logo.ico'))
            # ... or something went wrong
            except:
                raise Exception('croco_logo.ico not found')
        ## end setting the icon

        icon = wx.Icon(iconFile, wx.BITMAP_TYPE_ICO)
        self.SetIcon(icon)

        # generate a variable containing the parent to allow pushing back
        # the user-provided input
        self.parent = parent

        # set starting path for folder selection GUI
        self.currentPath = os.getcwd()

        # in these variables the options from CroCoMainWindow will be written
        self.InputOptionsToAsk = ''
        self.OutputOptionsToAsk = ''

    def updateOptions(self):

        def CreateOptionsContainer(self, option, listofOptions):
            """
            Args:
                option: tuple with (label, type of input)
                listOfOptions: a list object to append to user-input to
            Returns:
                optionContainer: wx.BoxSizer object
            """

            optionLabel = wx.StaticText(self.panel,
                                        wx.ID_ANY,
                                        label=option[0],
                                        style=wx.ALIGN_CENTER)
            optionLabel.Bind(wx.EVT_HELP,
                             lambda evt: self.Info(option[2], caption='Help'))

            if option[1] == 'file':
                optionButton = wx.Button(self.panel, label='Load file', name=option[0])
                # Bind using a lambda function to be able to pass multiple args
                # to the called function
                # see: https://wiki.wxpython.org/Passing%20Arguments%20to%20Callbacks
                optionButton.Bind(wx.EVT_BUTTON,
                                  lambda evt, appendTo=listofOptions, label=option[0]: self.onOpenFile(evt, appendTo, label))
                optionButton.Bind(wx.EVT_HELP,
                                  lambda evt: self.Info(option[2], caption='Help'))

            elif option[1] == 'dir':
                optionButton = wx.Button(self.panel, label='Open directory', name=option[0])
                optionButton.Bind(wx.EVT_BUTTON,
                                  lambda evt, appendTo=listofOptions, label=option[0]: self.onOpenDir(evt, appendTo, label))
                optionButton.Bind(wx.EVT_HELP,
                                  lambda evt: self.Info(option[2], caption='Help'))
            elif option[1] == 'input':
                # The 4th (optional) element of the options tuple is the multi
                # purpose field that may contain the default value for the
                # input field
                if len(option) > 3:
                    defValue = option[3]
                else:
                    defValue = ''
                optionButton = wx.TextCtrl(self.panel, value=defValue, name=option[0])
                self.textCtrls[option[0]] = optionButton
                optionButton.Bind(wx.EVT_HELP,
                                  lambda evt: self.Info(option[2], caption='Help'))
            elif option[1] == 'check':
                optionButton = wx.CheckBox(self.panel, name=option[0])
                self.textCtrls[option[0]] = optionButton
                optionButton.Bind(wx.EVT_HELP,
                                  lambda evt: self.Info(option[2], caption='Help'))
            else:
                raise Exception('[createOptionsContainer] Wrong options format for option: {}'\
                                    .format(option[0]))

            optionContainer = wx.BoxSizer(wx.HORIZONTAL)
            optionContainer.Add(optionLabel, 0, wx.ALL|wx.EXPAND, 5)
            optionContainer.Add(optionButton, 0, wx.ALL|wx.EXPAND, 5)

            return optionContainer

        # Dicts mapping the label names to the user-provided input
        self.inputOptionsToUserInput = {}
        self.outputOptionsToUserInput = {}

        # dict mapping text labels to input text fields
        self.textCtrls = {}

        # lists containign the label, button sizers for passing to higher
        # level sizer
        InputOptionContainers = []
        OutputOptionContainers = []

        for option in self.InputOptionsToAsk:

            optionContainer = CreateOptionsContainer(self,
                                                     option,
                                                     self.inputOptionsToUserInput)
            InputOptionContainers.append((optionContainer, 0, wx.ALL|wx.EXPAND, 5))

        for option in self.OutputOptionsToAsk:
            optionContainer = CreateOptionsContainer(self,
                                                     option,
                                                     self.outputOptionsToUserInput)
            OutputOptionContainers.append((optionContainer, 0, wx.ALL|wx.EXPAND, 5))

        ### createWidgets

        inputSizer = wx.BoxSizer(wx.VERTICAL)
        inputSizer.AddMany(InputOptionContainers)

        outputSizer = wx.BoxSizer(wx.VERTICAL)
        outputSizer.AddMany(OutputOptionContainers)

        optionSizer = wx.BoxSizer(wx.HORIZONTAL)
        optionSizer.Add(inputSizer, 0, wx.ALL|wx.EXPAND, 5)
        optionSizer.Add(outputSizer, 0, wx.ALL|wx.EXPAND, 5)

        ## Controls
        okayButton = wx.Button(self.panel, label='Okay')
        closeButton = wx.Button(self.panel, label='Close')
        helpButton = wx.ContextHelpButton(self.panel)

        controlSizer = wx.BoxSizer(wx.HORIZONTAL)
        controlSizer.Add(okayButton, 0, wx.EXPAND|wx.LEFT, 5)
        controlSizer.Add(closeButton, 0, wx.EXPAND|wx.RIGHT, 5)
        controlSizer.Add(helpButton, 0, wx.LEFT | wx.ALL | wx.EXPAND, 5)

        ## topSizer
        topSizer = wx.BoxSizer(wx.VERTICAL)
        topSizer.Add(optionSizer, 0, wx.ALL|wx.EXPAND, 5)
        topSizer.Add(wx.StaticLine(self.panel), 0, wx.ALL|wx.EXPAND, 5)
        topSizer.Add(controlSizer, 0, wx.ALL|wx.EXPAND, 5)

        # assign the top sizer to the panel i.e. main layout instance
        self.panel.SetSizer(topSizer)
        topSizer.Fit(self)

        closeButton.Bind(wx.EVT_BUTTON, self.onCancel)
        okayButton.Bind(wx.EVT_BUTTON, self.onOkay)

    def onCancel(self, event):
        print('[onCancel] Closing options window')
        self.Close()

    def onOkay(self, event):

        inputArgs = []
        outputArgs = []

        # label, type, help, optional
        for option in self.OutputOptionsToAsk:
            # all inputs done by extra dialog widgets (file and paths) are
            # stored upon closing
            # the dialogs. Text and checkbox inputs have to be evaaluated
            # when the okay-Button has been pressed
            if (option[1] == 'input') or (option[1] == 'check'):
                userInput = self.textCtrls[option[0]].GetValue()
                self.outputOptionsToUserInput[option[0]] = userInput
            userInput = self.outputOptionsToUserInput[option[0]]
            outputArgs.append(userInput)

        for option in self.InputOptionsToAsk:
            if (option[1] == 'input') or (option[1] == 'check'):
                userInput = self.textCtrls[option[0]].GetValue()
                self.inputOptionsToUserInput[option[0]] = userInput
            userInput = self.inputOptionsToUserInput[option[0]]
            inputArgs.append(userInput)

        print('[onOkay] Options for Input')
        for key in self.inputOptionsToUserInput:
            print('\t{}: {}'.format(key, self.inputOptionsToUserInput[key]))

        print('[onOkay] Options for Output')
        for key in self.outputOptionsToUserInput:
            print('\t{}: {}'.format(key, self.outputOptionsToUserInput[key]))

        self.parent.inputArgs = inputArgs
        self.parent.outputArgs = outputArgs

        self.parent.onRun(event)

        # Closing options window after passing options to main window
        self.Close()

    def onOpenFile(self, event, dictToAppend, label):
       dlg = wx.FileDialog(self,
                           message="Choose a file for input",
                           defaultDir=self.currentPath,
                           defaultFile="",
                           wildcard="*.*",
                           style=wx.FD_MULTIPLE)
       if dlg.ShowModal() == wx.ID_OK:
            self.currentPath = dictToAppend[label] = dlg.GetPath()
            print('[onOpenFile] Loaded {}'.format(self.currentPath))
       dlg.Destroy()

    def onOpenDir(self, event, dictToAppend, label):
        dlg = wx.DirDialog(self,
                           message="Choose one directory for Input:",
                           defaultPath=self.currentPath,
                           style=wx.DD_DEFAULT_STYLE | wx.DD_DIR_MUST_EXIST)

        if dlg.ShowModal() == wx.ID_OK:
            self.currentPath = dictToAppend[label] = dlg.GetPath()
            print('[onOpenDir] Loaded {}'.format(self.currentPath))
        dlg.Destroy()

    def Info(self, message, caption = 'CroCo'):
        dlg = wx.MessageDialog(self, message, caption, wx.OK | wx.ICON_INFORMATION)
        dlg.ShowModal()
        dlg.Destroy()


if __name__ == '__main__':
    # When this module is run (not imported) then create the app, the
    # frame, show it, and start the event loop.
    app = wx.App()
    frm = CroCoMainFrame().Show()
    app.MainLoop()