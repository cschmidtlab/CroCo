# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'qt\spectrumOptions.ui'
#
# Created by: PyQt5 UI code generator 5.6
#
# WARNING! All changes made in this file will be lost!

from PyQt5 import QtCore, QtGui, QtWidgets

class Ui_SpectrumAssignmentOptions(object):
    def setupUi(self, SpectrumAssignmentOptions):
        SpectrumAssignmentOptions.setObjectName("SpectrumAssignmentOptions")
        SpectrumAssignmentOptions.resize(362, 436)
        self.centralwidget = QtWidgets.QWidget(SpectrumAssignmentOptions)
        self.centralwidget.setObjectName("centralwidget")
        self.gridLayout = QtWidgets.QGridLayout(self.centralwidget)
        self.gridLayout.setObjectName("gridLayout")
        self.pushButton = QtWidgets.QPushButton(self.centralwidget)
        self.pushButton.setObjectName("pushButton")
        self.gridLayout.addWidget(self.pushButton, 3, 0, 1, 1)
        self.pushButton_2 = QtWidgets.QPushButton(self.centralwidget)
        self.pushButton_2.setObjectName("pushButton_2")
        self.gridLayout.addWidget(self.pushButton_2, 3, 1, 1, 1)
        self.groupBox_2 = QtWidgets.QGroupBox(self.centralwidget)
        self.groupBox_2.setObjectName("groupBox_2")
        self.gridLayout_2 = QtWidgets.QGridLayout(self.groupBox_2)
        self.gridLayout_2.setObjectName("gridLayout_2")
        self.label = QtWidgets.QLabel(self.groupBox_2)
        self.label.setObjectName("label")
        self.gridLayout_2.addWidget(self.label, 1, 0, 1, 1)
        self.gridLayout.addWidget(self.groupBox_2, 1, 0, 1, 2)
        self.groupBox_3 = QtWidgets.QGroupBox(self.centralwidget)
        self.groupBox_3.setObjectName("groupBox_3")
        self.formLayout_2 = QtWidgets.QFormLayout(self.groupBox_3)
        self.formLayout_2.setObjectName("formLayout_2")
        self.label_5 = QtWidgets.QLabel(self.groupBox_3)
        self.label_5.setObjectName("label_5")
        self.formLayout_2.setWidget(0, QtWidgets.QFormLayout.LabelRole, self.label_5)
        self.gridLayout.addWidget(self.groupBox_3, 2, 0, 1, 2)
        self.groupBox = QtWidgets.QGroupBox(self.centralwidget)
        self.groupBox.setObjectName("groupBox")
        self.formLayout = QtWidgets.QFormLayout(self.groupBox)
        self.formLayout.setObjectName("formLayout")
        self.label_2 = QtWidgets.QLabel(self.groupBox)
        self.label_2.setObjectName("label_2")
        self.formLayout.setWidget(0, QtWidgets.QFormLayout.LabelRole, self.label_2)
        self.mass_mono = QtWidgets.QRadioButton(self.groupBox)
        self.mass_mono.setChecked(True)
        self.mass_mono.setObjectName("mass_mono")
        self.options_ionmass = QtWidgets.QButtonGroup(SpectrumAssignmentOptions)
        self.options_ionmass.setObjectName("options_ionmass")
        self.options_ionmass.addButton(self.mass_mono)
        self.formLayout.setWidget(1, QtWidgets.QFormLayout.LabelRole, self.mass_mono)
        self.mass_av = QtWidgets.QRadioButton(self.groupBox)
        self.mass_av.setObjectName("mass_av")
        self.options_ionmass.addButton(self.mass_av)
        self.formLayout.setWidget(1, QtWidgets.QFormLayout.FieldRole, self.mass_av)
        self.label_4 = QtWidgets.QLabel(self.groupBox)
        self.label_4.setObjectName("label_4")
        self.formLayout.setWidget(2, QtWidgets.QFormLayout.LabelRole, self.label_4)
        self.options_aion = QtWidgets.QCheckBox(self.groupBox)
        self.options_aion.setChecked(True)
        self.options_aion.setObjectName("options_aion")
        self.formLayout.setWidget(3, QtWidgets.QFormLayout.LabelRole, self.options_aion)
        self.options_xion = QtWidgets.QCheckBox(self.groupBox)
        self.options_xion.setObjectName("options_xion")
        self.formLayout.setWidget(3, QtWidgets.QFormLayout.FieldRole, self.options_xion)
        self.options_yion = QtWidgets.QCheckBox(self.groupBox)
        self.options_yion.setChecked(True)
        self.options_yion.setObjectName("options_yion")
        self.formLayout.setWidget(4, QtWidgets.QFormLayout.FieldRole, self.options_yion)
        self.options_cion = QtWidgets.QCheckBox(self.groupBox)
        self.options_cion.setObjectName("options_cion")
        self.formLayout.setWidget(5, QtWidgets.QFormLayout.LabelRole, self.options_cion)
        self.options_zion = QtWidgets.QCheckBox(self.groupBox)
        self.options_zion.setObjectName("options_zion")
        self.formLayout.setWidget(5, QtWidgets.QFormLayout.FieldRole, self.options_zion)
        self.label_3 = QtWidgets.QLabel(self.groupBox)
        self.label_3.setObjectName("label_3")
        self.formLayout.setWidget(6, QtWidgets.QFormLayout.LabelRole, self.label_3)
        self.lineEdit_2 = QtWidgets.QLineEdit(self.groupBox)
        self.lineEdit_2.setObjectName("lineEdit_2")
        self.formLayout.setWidget(7, QtWidgets.QFormLayout.LabelRole, self.lineEdit_2)
        self.options_bion = QtWidgets.QCheckBox(self.groupBox)
        self.options_bion.setChecked(True)
        self.options_bion.setObjectName("options_bion")
        self.formLayout.setWidget(4, QtWidgets.QFormLayout.LabelRole, self.options_bion)
        self.gridLayout.addWidget(self.groupBox, 0, 0, 1, 2)
        SpectrumAssignmentOptions.setCentralWidget(self.centralwidget)
        self.menubar = QtWidgets.QMenuBar(SpectrumAssignmentOptions)
        self.menubar.setGeometry(QtCore.QRect(0, 0, 362, 19))
        self.menubar.setObjectName("menubar")
        SpectrumAssignmentOptions.setMenuBar(self.menubar)

        self.retranslateUi(SpectrumAssignmentOptions)
        QtCore.QMetaObject.connectSlotsByName(SpectrumAssignmentOptions)

    def retranslateUi(self, SpectrumAssignmentOptions):
        _translate = QtCore.QCoreApplication.translate
        SpectrumAssignmentOptions.setWindowTitle(_translate("SpectrumAssignmentOptions", "Assignment Options"))
        self.pushButton.setText(_translate("SpectrumAssignmentOptions", "Save"))
        self.pushButton_2.setText(_translate("SpectrumAssignmentOptions", "Cancel"))
        self.groupBox_2.setTitle(_translate("SpectrumAssignmentOptions", "Cross-linker"))
        self.label.setText(_translate("SpectrumAssignmentOptions", "Not implemented! :D"))
        self.groupBox_3.setTitle(_translate("SpectrumAssignmentOptions", "Modifications"))
        self.label_5.setText(_translate("SpectrumAssignmentOptions", "Not implemented, too!"))
        self.groupBox.setTitle(_translate("SpectrumAssignmentOptions", "General"))
        self.label_2.setText(_translate("SpectrumAssignmentOptions", "Mass of the ions"))
        self.mass_mono.setText(_translate("SpectrumAssignmentOptions", "monisotopic"))
        self.mass_av.setText(_translate("SpectrumAssignmentOptions", "average"))
        self.label_4.setText(_translate("SpectrumAssignmentOptions", "Ion types"))
        self.options_aion.setText(_translate("SpectrumAssignmentOptions", "a"))
        self.options_xion.setText(_translate("SpectrumAssignmentOptions", "x"))
        self.options_yion.setText(_translate("SpectrumAssignmentOptions", "y"))
        self.options_cion.setText(_translate("SpectrumAssignmentOptions", "c"))
        self.options_zion.setText(_translate("SpectrumAssignmentOptions", "z"))
        self.label_3.setText(_translate("SpectrumAssignmentOptions", "Max m/z to consider (depends on quad)"))
        self.lineEdit_2.setText(_translate("SpectrumAssignmentOptions", "2000"))
        self.options_bion.setText(_translate("SpectrumAssignmentOptions", "b"))


if __name__ == "__main__":
    import sys
    app = QtWidgets.QApplication(sys.argv)
    SpectrumAssignmentOptions = QtWidgets.QMainWindow()
    ui = Ui_SpectrumAssignmentOptions()
    ui.setupUi(SpectrumAssignmentOptions)
    SpectrumAssignmentOptions.show()
    sys.exit(app.exec_())

