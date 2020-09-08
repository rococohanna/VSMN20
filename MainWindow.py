# -*- coding: utf-8 -*-
"""
Created on Tue May  7 09:12:27 2019

@author: Ludvig Hassbring & Hanna Nilsson
"""

# -*- coding: utf-8 -*-

import sys

from PyQt5.QtCore import pyqtSlot, pyqtSignal, QThread
from PyQt5.QtWidgets import QApplication, QDialog, QWidget, QMainWindow, QFileDialog,QMessageBox,QPushButton
from PyQt5.uic import loadUi
from PyQt5 import QtGui, QtCore

import plainstress as ps
import json
import calfem.ui as cfui
import numpy as np




class SolverThread(QThread):
    """Klass för att hantera beräkning i bakgrunden"""

    def __init__(self, solver, paramStudy):
        """Klasskonstruktor"""
        QThread.__init__(self)
        self.solver = solver
        self.paramStudy = paramStudy

    def __del__(self):
        self.wait()

    def run(self):
       
        if(self.paramStudy):
            self.solver.executeParamStudy()
        
        else:
        #exe = ps.Solver(self.inputData,self.outputData)
            self.solver.execute()
        


class MainWindow(QMainWindow):
    """MainWindow-klass som hanterar vårt huvudfönster"""

    def __init__(self):
        """Constructor"""
        super(QMainWindow, self).__init__()
        self.calcDone = "false" #Flagga för att se till att det inte skiter sig 
        # --- Lagra en referens till applikationsinstansen i klassen

        self.app = app
        self.filename =""
        # --- Läs in gränssnitt från fil

        self.ui = loadUi('mainwindow.ui', self)

        # --- Se till att visa fönstret
        self.ui.show()
        self.ui.raise_()
        
        self.initModel()
        self.updateControls()
        self.calcDone  = False
        
        # --- Koppla kontroller till händelsemetoder

        self.ui.actionNew.triggered.connect(self.onActionNew)
        self.ui.actionOpen.triggered.connect(self.onActionOpen)
        self.ui.actionSave.triggered.connect(self.onActionSave)
        self.ui.actionExit.triggered.connect(self.onActionExit)
        self.ui.actionSave_as.triggered.connect(self.onActionSave_as)
        self.ui.actionExecute.triggered.connect(self.onActionExecute)
        
        
        
        
        # Skapar Vis-objekt första gången 
        
        self.vis = ps.Visualisation(self.inputData,self.outputData,self.calcDone)
        
        #Flagga finns i PlainStress
        self.ui.showElementValuesButton.clicked.connect(self.showElementValues)
        self.ui.showNodalValuesButton.clicked.connect(self.showNodalValues)
        self.ui.showGeometryButton.clicked.connect(self.showGeometry)
        self.ui.showMeshButton.clicked.connect(self.showMesh)
        self.ui.paramButton.clicked.connect(self.onExecuteParamStudy)
        
        
        
        
      



    def onActionNew(self):
        """Skapa en ny modell"""
        self.filename = ""
        self.ui.plainTextEdit.setPlainText("")
        
        
    
    def onActionExit(self):
        """Avslutar Programmet"""
        choice = QMessageBox.question(self, 'Exit',
                                            "Would you like to exit the program?",
                                            QMessageBox.Yes | QMessageBox.No)
        if choice == QMessageBox.Yes:
            print("Shutting Down..")
            sys.exit()
        else:
            pass
        
    def onActionSave_as(self):
        self.filename, _  = QFileDialog.getSaveFileName(self.ui, 
                "Spara modell", "", "Modell filer (*.json)")
        """Skapa en ny modell"""
        inputData = {}
        inputData["version"] = self.savedata.version
        inputData["t"] = self.savedata.t
        inputData["ep"] = [self.savedata.ep[0],self.savedata.ep[1]]
        inputData["w"] = self.savedata.w
        inputData["a"] = self.savedata.t
        inputData["b"] = self.savedata.b
        inputData["h"] = self.savedata.h
        inputData["v"] = self.savedata.v
        inputData["E"] = self.savedata.E
        ofile = open(self.filename, "w")
        json.dump(inputData, ofile, sort_keys = True, indent = 4)
        ofile.close()
    
        
    def onActionOpen(self):
         """Öppna in indata fil"""

         self.filename, _ = QFileDialog.getOpenFileName(self.ui, "Öppna modell", "", "Modell filer (*.json *.jpg *.bmp)")
        

         if self.filename!="":
             ifile = open(self.filename, "r")
             inputData = json.load(ifile)
             ifile.close()
    
             self.inputData.version = inputData["h"]
             self.inputData.bcs = np.asarray(inputData["b"])
             self.inputData.t = inputData["t"]
             self.inputData.ep = np.asarray(inputData["ep"])
             self.inputData.dof = np.asarray(inputData["a"])
             self.inputData.w = inputData["w"]
             self.inputData.E =inputData["E"]
             self.inputData.v = inputData["v"]
             self.inputData.Elementsize = inputData["Elementsize"]
             self.updateControls()
        
    def onActionSave(self):
        """Spara modell"""
    
        self.updateModel()
        
        self.savedata = ps.InputData()
        if self.filename == "":
            self.filename, _  = QFileDialog.getSaveFileName(self.ui, 
                "Spara modell", "", "Modell filer (*.json)")
    
        if self.filename!="":
            inputData = {}
            inputData["version"] = self.savedata.version
            inputData["t"] = self.savedata.t
            inputData["ep"] = [self.savedata.ep[0],self.savedata.ep[1]]
            inputData["w"] = self.savedata.w
            inputData["a"] = self.savedata.t
            inputData["b"] = self.savedata.b
            inputData["h"] = self.savedata.h
            inputData["Elementsize"] = self.savedata.Elementsize
            inputData["v"] = self.savedata.v
            inputData["E"] = self.savedata.E
            ofile = open(self.filename, "w")
            json.dump(inputData, ofile, sort_keys = True, indent = 4)
            ofile.close()
    def onActionExecute(self):
        """Kör beräkningen"""

        # --- Avaktivera gränssnitt under beräkningen.
        
        
        self.ui.setEnabled(False)

        # --- Uppdatera värden från kontroller
   
        self.updateModel()
        
        # --- Återställ värden
        
        
        ps.Visualisation.closeAll(self)

        # --- Skapa en lösare

        self.solver = ps.Solver(self.inputData, self.outputData)

        # --- Starta en tråd för att köra beräkningen, så att 
        #     gränssnittet inte fryser.

        self.solverThread = SolverThread(self.solver,paramStudy=False)     
        self.solverThread.finished.connect(self.onSolverFinished)  
        self.solverThread.start()      
        
        print(self.inputData.Elementsize)
        self.vis.closeAll()
        
        
    def onSolverFinished(self):
        """Anropas när beräkningstråden avslutas"""

        # --- Aktivera gränssnitt igen
        
        

        self.ui.setEnabled(True)

        # --- Generera resulatrappor
        self.calcDone = True
        self.vis = ps.Visualisation(self.inputData,self.outputData,self.calcDone)
      
        self.ui.plainTextEdit.setPlainText(str(ps.Report(self.inputData,self.outputData)))
    def showGeometry(self):
        self.vis.showGeometry()
    def showElementValues(self):
        self.vis.showElementValues()
    def showNodalValues(self):
        self.vis.showNodalValues()
    def showMesh(self):
        self.vis.showMesh()
    def initModel(self):
        """Initierar värden på modellen"""
     
       
        self.inputData = ps.InputData()
        self.outputData= ps.OutputData()
        
        
    def updateControls(self):
        """Fyll kontrollerna med värden från modellen"""
        
        self.ui.wEdit.setText(str(self.inputData.w))
        self.ui.hEdit.setText(str(self.inputData.h))
        self.ui.tEdit.setText(str(self.inputData.t))
        self.ui.dEdit.setText(str(self.inputData.a))
        
        self.ui.vEdit.setText(str(self.inputData.v))
        self.ui.eEdit.setText(str(self.inputData.E))
        
        #self.Elementsize = 0.01
        #self.ui.horizontalSlider.setValue(int(self.inputData.Elementsize))*1000
      
        
    
    def updateModel(self):
        """Hämta värden från kontroller och uppdatera modellen"""

        self.inputData.w = float(self.ui.wEdit.text())
        self.inputData.t = float(self.ui.tEdit.text())
        self.inputData.a = float(self.ui.dEdit.text())
        self.inputData.h = float(self.ui.hEdit.text())
        
        self.inputData.E = float(self.ui.eEdit.text())
        self.inputData.v = float(self.ui.vEdit.text())
        self.inputData.Elementsize = float(self.ui.horizontalSlider.value())/1000
        
        
    def onExecuteParamStudy(self):
        """Exekvera parameterstudie"""
    
        # --- Hämta värden från grafiskt gränssnitt.
        
    
        self.inputData.paramD = self.ui.paramVaryDRadio.isChecked() #dubbelkolla namn
        self.inputData.paramT = self.ui.paramVaryTRadio.isChecked()
    
        if self.inputData.paramD:
            self.inputData.dStart = float(self.ui.dEdit.text())
            self.inputData.dEnd = float(self.ui.dEndEdit.text())
        elif self.inputData.paramT:
            self.inputData.tStart = float(self.ui.tEdit.text())
            self.inputData.tEnd = float(self.ui.tEndEdit.text())
    
        self.inputData.paramFilename = "paramStudy"
        self.inputData.paramSteps = int(self.ui.paramStep.value())
        
            # --- Starta en tråd för att köra beräkningen, så att 
            #     gränssnittet inte fryser.
        self.solver = ps.Solver(self.inputData, self.outputData)
        self.solverThread = SolverThread(self.solver, paramStudy = True)        
        self.solverThread.finished.connect(self.onSolverFinished)        
        self.solverThread.start()
        
    
    
if __name__ == '__main__':

    # --- Skapa applikationsinstans

    app = QApplication(sys.argv)

    # --- Skapa och visa huvudfönster

    widget = MainWindow()
    widget.show()

    # --- Starta händelseloopen
    
    sys.exit(app.exec_())