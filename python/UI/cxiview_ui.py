# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'cxiview.ui'
#
# Created by: PyQt5 UI code generator 5.9.2
#
# WARNING! All changes made in this file will be lost!

from PyQt5 import QtCore, QtGui, QtWidgets

class Ui_MainWindow(object):
    def setupUi(self, MainWindow):
        MainWindow.setObjectName("MainWindow")
        MainWindow.resize(1400, 1005)
        self.centralwidget = QtWidgets.QWidget(MainWindow)
        self.centralwidget.setObjectName("centralwidget")
        self.gridLayout = QtWidgets.QGridLayout(self.centralwidget)
        self.gridLayout.setObjectName("gridLayout")
        self.verticalLayout = QtWidgets.QVBoxLayout()
        self.verticalLayout.setObjectName("verticalLayout")
        self.horizontalLayout4 = QtWidgets.QHBoxLayout()
        self.horizontalLayout4.setObjectName("horizontalLayout4")
        self.verticalLayout.addLayout(self.horizontalLayout4)
        self.horizontalLayout = QtWidgets.QHBoxLayout()
        self.horizontalLayout.setSizeConstraint(QtWidgets.QLayout.SetFixedSize)
        self.horizontalLayout.setObjectName("horizontalLayout")
        self.refreshfilesPushButton = QtWidgets.QPushButton(self.centralwidget)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Fixed, QtWidgets.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.refreshfilesPushButton.sizePolicy().hasHeightForWidth())
        self.refreshfilesPushButton.setSizePolicy(sizePolicy)
        self.refreshfilesPushButton.setMaximumSize(QtCore.QSize(100, 32))
        self.refreshfilesPushButton.setObjectName("refreshfilesPushButton")
        self.horizontalLayout.addWidget(self.refreshfilesPushButton)
        self.previousPushButton = QtWidgets.QPushButton(self.centralwidget)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Fixed, QtWidgets.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.previousPushButton.sizePolicy().hasHeightForWidth())
        self.previousPushButton.setSizePolicy(sizePolicy)
        self.previousPushButton.setMaximumSize(QtCore.QSize(100, 16777215))
        self.previousPushButton.setObjectName("previousPushButton")
        self.horizontalLayout.addWidget(self.previousPushButton)
        self.nextPushButton = QtWidgets.QPushButton(self.centralwidget)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Fixed, QtWidgets.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.nextPushButton.sizePolicy().hasHeightForWidth())
        self.nextPushButton.setSizePolicy(sizePolicy)
        self.nextPushButton.setMaximumSize(QtCore.QSize(100, 16777215))
        self.nextPushButton.setObjectName("nextPushButton")
        self.horizontalLayout.addWidget(self.nextPushButton)
        self.playPushButton = QtWidgets.QPushButton(self.centralwidget)
        self.playPushButton.setEnabled(True)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Fixed, QtWidgets.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.playPushButton.sizePolicy().hasHeightForWidth())
        self.playPushButton.setSizePolicy(sizePolicy)
        self.playPushButton.setMaximumSize(QtCore.QSize(100, 16777215))
        self.playPushButton.setObjectName("playPushButton")
        self.horizontalLayout.addWidget(self.playPushButton)
        self.randomPushButton = QtWidgets.QPushButton(self.centralwidget)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Fixed, QtWidgets.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.randomPushButton.sizePolicy().hasHeightForWidth())
        self.randomPushButton.setSizePolicy(sizePolicy)
        self.randomPushButton.setMaximumSize(QtCore.QSize(100, 16777215))
        self.randomPushButton.setObjectName("randomPushButton")
        self.horizontalLayout.addWidget(self.randomPushButton)
        self.shufflePushButton = QtWidgets.QPushButton(self.centralwidget)
        self.shufflePushButton.setEnabled(True)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Fixed, QtWidgets.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.shufflePushButton.sizePolicy().hasHeightForWidth())
        self.shufflePushButton.setSizePolicy(sizePolicy)
        self.shufflePushButton.setMaximumSize(QtCore.QSize(100, 16777215))
        self.shufflePushButton.setObjectName("shufflePushButton")
        self.horizontalLayout.addWidget(self.shufflePushButton)
        self.jumpToLabel = QtWidgets.QLabel(self.centralwidget)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Fixed, QtWidgets.QSizePolicy.Preferred)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.jumpToLabel.sizePolicy().hasHeightForWidth())
        self.jumpToLabel.setSizePolicy(sizePolicy)
        self.jumpToLabel.setMaximumSize(QtCore.QSize(100, 16777215))
        self.jumpToLabel.setAlignment(QtCore.Qt.AlignRight|QtCore.Qt.AlignTrailing|QtCore.Qt.AlignVCenter)
        self.jumpToLabel.setObjectName("jumpToLabel")
        self.horizontalLayout.addWidget(self.jumpToLabel)
        self.jumpToLineEdit = QtWidgets.QLineEdit(self.centralwidget)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Fixed, QtWidgets.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.jumpToLineEdit.sizePolicy().hasHeightForWidth())
        self.jumpToLineEdit.setSizePolicy(sizePolicy)
        self.jumpToLineEdit.setMaximumSize(QtCore.QSize(100, 16777215))
        self.jumpToLineEdit.setObjectName("jumpToLineEdit")
        self.horizontalLayout.addWidget(self.jumpToLineEdit)
        spacerItem = QtWidgets.QSpacerItem(40, 20, QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Minimum)
        self.horizontalLayout.addItem(spacerItem)
        self.verticalLayout.addLayout(self.horizontalLayout)
        self.imageView = ImageView(self.centralwidget)
        self.imageView.setEnabled(True)
        self.imageView.setObjectName("imageView")
        self.verticalLayout.addWidget(self.imageView)
        self.horizontalLayout2 = QtWidgets.QHBoxLayout()
        self.horizontalLayout2.setObjectName("horizontalLayout2")
        self.horizontalLayout3 = QtWidgets.QHBoxLayout()
        self.horizontalLayout3.setObjectName("horizontalLayout3")
        self.statusBar = QtWidgets.QLabel(self.centralwidget)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.Preferred)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.statusBar.sizePolicy().hasHeightForWidth())
        self.statusBar.setSizePolicy(sizePolicy)
        self.statusBar.setMinimumSize(QtCore.QSize(500, 0))
        self.statusBar.setObjectName("statusBar")
        self.horizontalLayout3.addWidget(self.statusBar)
        spacerItem1 = QtWidgets.QSpacerItem(40, 20, QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Minimum)
        self.horizontalLayout3.addItem(spacerItem1)
        self.foundPeaksCheckBox = QtWidgets.QCheckBox(self.centralwidget)
        self.foundPeaksCheckBox.setMaximumSize(QtCore.QSize(150, 16777215))
        self.foundPeaksCheckBox.setObjectName("foundPeaksCheckBox")
        self.horizontalLayout3.addWidget(self.foundPeaksCheckBox)
        self.predictedPeaksCheckBox = QtWidgets.QCheckBox(self.centralwidget)
        self.predictedPeaksCheckBox.setEnabled(False)
        self.predictedPeaksCheckBox.setMaximumSize(QtCore.QSize(150, 16777215))
        self.predictedPeaksCheckBox.setObjectName("predictedPeaksCheckBox")
        self.horizontalLayout3.addWidget(self.predictedPeaksCheckBox)
        self.masksCheckBox = QtWidgets.QCheckBox(self.centralwidget)
        self.masksCheckBox.setEnabled(True)
        self.masksCheckBox.setMaximumSize(QtCore.QSize(150, 16777215))
        self.masksCheckBox.setObjectName("masksCheckBox")
        self.horizontalLayout3.addWidget(self.masksCheckBox)
        self.resolutionCheckBox = QtWidgets.QCheckBox(self.centralwidget)
        self.resolutionCheckBox.setEnabled(True)
        self.resolutionCheckBox.setMaximumSize(QtCore.QSize(150, 16777215))
        self.resolutionCheckBox.setObjectName("resolutionCheckBox")
        self.horizontalLayout3.addWidget(self.resolutionCheckBox)
        self.horizontalLayout2.addLayout(self.horizontalLayout3)
        self.verticalLayout.addLayout(self.horizontalLayout2)
        self.gridLayout.addLayout(self.verticalLayout, 1, 0, 1, 1)
        MainWindow.setCentralWidget(self.centralwidget)
        self.menuBar = QtWidgets.QMenuBar(MainWindow)
        self.menuBar.setGeometry(QtCore.QRect(0, 0, 1400, 22))
        self.menuBar.setObjectName("menuBar")
        self.menuFile = QtWidgets.QMenu(self.menuBar)
        self.menuFile.setObjectName("menuFile")
        self.menuColours = QtWidgets.QMenu(self.menuBar)
        self.menuColours.setObjectName("menuColours")
        self.menuAnalysis = QtWidgets.QMenu(self.menuBar)
        self.menuAnalysis.setObjectName("menuAnalysis")
        self.menuCXI = QtWidgets.QMenu(self.menuBar)
        self.menuCXI.setObjectName("menuCXI")
        self.menuCrystals = QtWidgets.QMenu(self.menuBar)
        self.menuCrystals.setObjectName("menuCrystals")
        self.menuView = QtWidgets.QMenu(self.menuBar)
        self.menuView.setObjectName("menuView")
        MainWindow.setMenuBar(self.menuBar)
        self.actionSave_data = QtWidgets.QAction(MainWindow)
        self.actionSave_data.setEnabled(True)
        self.actionSave_data.setObjectName("actionSave_data")
        self.actionSave_image = QtWidgets.QAction(MainWindow)
        self.actionSave_image.setObjectName("actionSave_image")
        self.actionLoad_geometry = QtWidgets.QAction(MainWindow)
        self.actionLoad_geometry.setObjectName("actionLoad_geometry")
        self.actionRefresh_file_list = QtWidgets.QAction(MainWindow)
        self.actionRefresh_file_list.setObjectName("actionRefresh_file_list")
        self.actionQuit = QtWidgets.QAction(MainWindow)
        self.actionQuit.setObjectName("actionQuit")
        self.actionDefault_crystal_display_settings = QtWidgets.QAction(MainWindow)
        self.actionDefault_crystal_display_settings.setObjectName("actionDefault_crystal_display_settings")
        self.actionCircle_Cheetah_peaks = QtWidgets.QAction(MainWindow)
        self.actionCircle_Cheetah_peaks.setObjectName("actionCircle_Cheetah_peaks")
        self.actionDefault_particle_display_settings = QtWidgets.QAction(MainWindow)
        self.actionDefault_particle_display_settings.setObjectName("actionDefault_particle_display_settings")
        self.actionHistogram_clip = QtWidgets.QAction(MainWindow)
        self.actionHistogram_clip.setCheckable(True)
        self.actionHistogram_clip.setChecked(True)
        self.actionHistogram_clip.setObjectName("actionHistogram_clip")
        self.actionAuto_scale_levels = QtWidgets.QAction(MainWindow)
        self.actionAuto_scale_levels.setCheckable(True)
        self.actionAuto_scale_levels.setChecked(True)
        self.actionAuto_scale_levels.setObjectName("actionAuto_scale_levels")
        self.actionAutoscale = QtWidgets.QAction(MainWindow)
        self.actionAutoscale.setCheckable(True)
        self.actionAutoscale.setChecked(True)
        self.actionAutoscale.setObjectName("actionAutoscale")
        self.menu_view_photonconversion = QtWidgets.QAction(MainWindow)
        self.menu_view_photonconversion.setObjectName("menu_view_photonconversion")
        self.actionSave_data_assembled = QtWidgets.QAction(MainWindow)
        self.actionSave_data_assembled.setObjectName("actionSave_data_assembled")
        self.action_Imagefloorzero = QtWidgets.QAction(MainWindow)
        self.action_Imagefloorzero.setCheckable(True)
        self.action_Imagefloorzero.setChecked(True)
        self.action_Imagefloorzero.setObjectName("action_Imagefloorzero")
        self.menuFile.addAction(self.actionRefresh_file_list)
        self.menuFile.addAction(self.actionLoad_geometry)
        self.menuFile.addAction(self.actionSave_image)
        self.menuFile.addAction(self.actionSave_data)
        self.menuFile.addAction(self.actionSave_data_assembled)
        self.menuFile.addSeparator()
        self.menuFile.addAction(self.actionQuit)
        self.menuCXI.addAction(self.actionDefault_particle_display_settings)
        self.menuCrystals.addAction(self.actionDefault_crystal_display_settings)
        self.menuCrystals.addAction(self.actionCircle_Cheetah_peaks)
        self.menuView.addAction(self.actionAutoscale)
        self.menuView.addAction(self.actionAuto_scale_levels)
        self.menuView.addAction(self.actionHistogram_clip)
        self.menuView.addAction(self.action_Imagefloorzero)
        self.menuView.addAction(self.menu_view_photonconversion)
        self.menuBar.addAction(self.menuFile.menuAction())
        self.menuBar.addAction(self.menuColours.menuAction())
        self.menuBar.addAction(self.menuCrystals.menuAction())
        self.menuBar.addAction(self.menuCXI.menuAction())
        self.menuBar.addAction(self.menuAnalysis.menuAction())
        self.menuBar.addAction(self.menuView.menuAction())

        self.retranslateUi(MainWindow)
        QtCore.QMetaObject.connectSlotsByName(MainWindow)

    def retranslateUi(self, MainWindow):
        _translate = QtCore.QCoreApplication.translate
        MainWindow.setWindowTitle(_translate("MainWindow", "MainWindow"))
        self.refreshfilesPushButton.setText(_translate("MainWindow", "Refresh files"))
        self.previousPushButton.setText(_translate("MainWindow", "Previous"))
        self.nextPushButton.setText(_translate("MainWindow", "Next "))
        self.playPushButton.setText(_translate("MainWindow", "Play"))
        self.randomPushButton.setText(_translate("MainWindow", "Random"))
        self.shufflePushButton.setText(_translate("MainWindow", "Shuffle"))
        self.jumpToLabel.setText(_translate("MainWindow", "Frame"))
        self.statusBar.setText(_translate("MainWindow", "Last clicked pixel:     x: -     y: -     value: -     resolution: -"))
        self.foundPeaksCheckBox.setText(_translate("MainWindow", "Found peaks"))
        self.predictedPeaksCheckBox.setText(_translate("MainWindow", "Predicted peaks"))
        self.masksCheckBox.setText(_translate("MainWindow", "Pixel masks"))
        self.resolutionCheckBox.setText(_translate("MainWindow", "Resolution rings"))
        self.menuFile.setTitle(_translate("MainWindow", "File"))
        self.menuColours.setTitle(_translate("MainWindow", "Colours"))
        self.menuAnalysis.setTitle(_translate("MainWindow", "Analysis"))
        self.menuCXI.setTitle(_translate("MainWindow", "Particles"))
        self.menuCrystals.setTitle(_translate("MainWindow", "Crystals"))
        self.menuView.setTitle(_translate("MainWindow", "View"))
        self.actionSave_data.setText(_translate("MainWindow", "Save data (raw)"))
        self.actionSave_image.setText(_translate("MainWindow", "Save image"))
        self.actionLoad_geometry.setText(_translate("MainWindow", "Load geometry"))
        self.actionRefresh_file_list.setText(_translate("MainWindow", "Refresh file list"))
        self.actionQuit.setText(_translate("MainWindow", "Quit"))
        self.actionDefault_crystal_display_settings.setText(_translate("MainWindow", "Default crystal display settings"))
        self.actionCircle_Cheetah_peaks.setText(_translate("MainWindow", "Circle Cheetah peaks"))
        self.actionDefault_particle_display_settings.setText(_translate("MainWindow", "Default particle display settings"))
        self.actionHistogram_clip.setText(_translate("MainWindow", "Histogram clip"))
        self.actionAuto_scale_levels.setText(_translate("MainWindow", "Lock histogram scale"))
        self.actionAutoscale.setText(_translate("MainWindow", "Auto-scale image"))
        self.menu_view_photonconversion.setText(_translate("MainWindow", "Photon count conversion"))
        self.actionSave_data_assembled.setText(_translate("MainWindow", "Save data (assembled)"))
        self.action_Imagefloorzero.setText(_translate("MainWindow", "Image floor is zero"))

from pyqtgraph import ImageView
