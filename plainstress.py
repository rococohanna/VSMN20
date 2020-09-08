# -*- coding: utf-8 -*-
"""
Created on Thu Apr  4 13:28:54 2019

@author: Ludvig Hassbring & Hanna Nilsson 
"""



import numpy as np
import calfem.core as cfc
import calfem.geometry as cfg  # <-- Geometrirutiner
import calfem.mesh as cfm      # <-- Nätgenerering
import calfem.vis as cfv       # <-- Visualisering
import calfem.utils as cfu     # <-- Blandade rutiner
import pyvtk as vtk            # <-- Filexportering



import json

class InputData(object):
    """Klass för att definiera indata för vår modell."""
    def __init__(self):
       
        
        self.ptype = 1 #1=plain strain 
        self.h = 0.1
        self.t = 0.15
        self.w = 0.3
        self.a = 0.05
        self.b = 0.025
        self.E = 2.08e10 
        self.ep = [1,0.15]
        self.Elementsize = 0.01
        self.v = 0.2
        self.version = 1
       
        

    def geometry(self):
        """Skapa en geometri instans baserat på definierade parametrar"""

        # --- Skapa en geometri-instans för att lagra vår
        #     geometribeskrivning

        g = cfg.Geometry()

      
   
        
           # --- Enklare att använda referenser till self.xxx
        w = self.w
        h = self.h
        b = self.b
        a = self.a
        

        # --- Punkter i modellen skapas med point(...) metoden

        g.point([0, 0])
        g.point([(w-a)*0.5, 0])
        g.point([(w-a)*0.5, b])
        g.point([(w+a)*0.5, b])
        g.point([(w+a)*0.5, 0])
        g.point([w, 0])
        g.point([w, h])
        g.point([(w+a)*0.5, h])
        g.point([(w+a)*0.5, h-b])
        g.point([(w-a)*0.5, h-b])
        g.point([(w-a)*0.5, h])
        g.point([0, h])
        
       

        # --- Linker och kurvor skapas med spline(...) metoden


        q = 20
        wall = 30
        
        g.spline([0, 1])            
        g.spline([1, 2])           
        g.spline([2, 3])
        g.spline([3, 4])
        g.spline([4, 5])
        g.spline([5, 6],marker = q)
        g.spline([6, 7])
        g.spline([7, 8])
        g.spline([8, 9])
        g.spline([9,10])
        g.spline([10, 11])
        g.spline([11, 0],marker = wall)
       

        # --- Ytan på vilket nätet skall genereras definieras med 
        #     surface(...) metoden.

        g.surface([0,1,2,3,4,5,6,7,8,9,10,11])
        
        
        # --- Slutligen returnerar vi den skapade geometrin
        
       
        return g
        
    def save(self, filename):
        """Spara indata till fil."""

        inputData = {}
        inputData["version"] = self.version
        inputData["t"] = self.t
        inputData["ep"] = self.ep.tolist()
        inputData["w"] = self.w
        inputData["a"] = self.a
        inputData["b"] = self.b
        inputData["h"] = self.h
        inputData["v"] = self.v
        inputData["E"] = self.E
        

        ofile = open(filename, "w")
        json.dump(inputData, ofile, sort_keys = True, indent = 4)
        ofile.close()

    def load(self, filename):
        """Läs indata från fil."""

        ifile = open(filename, "r")
        inputData = json.load(ifile)
        ifile.close()

        self.version = inputData["h"]
        self.bcs = np.asarray(inputData["b"])
        self.t = inputData["t"]
        self.ep = np.asarray(inputData["ep"])
        self.dof = np.asarray(inputData["a"])
        self.w = inputData["w"]
        self.E =inputData["E"]
        self.v = inputData["v"]
        
          
        
        
class OutputData(object):
    """Klass för att lagra resultaten från beräkningen."""
    def __init__(self):
        self.asolve = None  #nodegenskap
        self.r = None #Reaktionskraft
        self.ed = None #Elementförskjutning
        self.qs = None #Elementflöden/spänningar
        self.qt = None #typ samma 
        self.vonMises = None #Maximala spänningar (Max av qs/qt?)
        self.edof = None #topology
        self.coords = None#coordinates
        self.elType = None#element type
        self.dofsPerNode = None#Degree of freedom per node
        self.geometry = None #geometry
        self.calcDone = False #Flagga för att se till att gränssnittet inte ritar utan execute 
        self.topt = None #Topologi
        self.stresses1 = None
        self.stresses2 = None 
        
class Solver(object):
    """Klass för att hantera lösningen av vår beräkningsmodell."""
    def __init__(self, inputData, outputData):
        self.inputData = inputData
        self.outputData = outputData
        
        
    def execute(self):

        # --- Överför modell variabler till lokala referenser

   
        ep = self.inputData.ep
        E = self.inputData.E
        v = self.inputData.v
        Elementsize = self.inputData.Elementsize
      
        
         # --- Anropa InputData för en geomtetribeskrivning

        geometry = self.inputData.geometry()
       
        # --- Nätgenerering

        elType = 3      # <-- Fyrnodselement flw2i4e 
        dofsPerNode= 2  

        meshGen = cfm.GmshMeshGenerator(geometry)
        
      
        meshGen.elSizeFactor = Elementsize     # <-- Anger max area för element
        meshGen.elType = elType
        meshGen.dofsPerNode = dofsPerNode
        meshGen.returnBoundaryElements = True

        coords, edof, dof, bdofs, elementmarkers, boundaryElements = meshGen.create()
        self.outputData.topo = meshGen.topo
        
        
        
        #Solver
        bc = np.array([],'i')
        bcVal = np.array([],'i')
       
      
    
        D = cfc.hooke(1,E,v)
        nDofs = np.size(dof)
        ex, ey = cfc.coordxtr(edof, coords, dof) #Coordinates 
        K = np.zeros([nDofs,nDofs])
        
        #Append Boundary Conds
        f = np.zeros([nDofs,1])
        bc, bcVal = cfu.applybc(bdofs, bc, bcVal, 30, 0.0, 0) 
        cfu.applyforce(bdofs, f, 20, 100e3, 1)
        
        
        qs_array = []
        qt_array = []
        
      
        
        for x,y,z in zip(ex,ey,edof):
            
            Ke= cfc.planqe(x,y,ep,D)
            
           
            cfc.assem(z, K, Ke)
            
           
       
        
        asolve,r = cfc.solveq(K,f,bc,bcVal)
       
        ed = cfc.extractEldisp(edof,asolve)
        for x,y,z in zip(ex,ey,ed):
            qs,qt = cfc.planqs(x,y,ep,D,z)
          
            qs_array.append(qs)
            qt_array.append(qt)
            
         
        vonMises = []
        stresses1 = []
        stresses2 = []
        # For each element:

        
          
        for i in range(edof.shape[0]): 

                # Determine element stresses and strains in the element.
    
            es, et = cfc.planqs(ex[i,:], ey[i,:], ep, D, ed[i,:]) 

            # Calc and append effective stress to list.    
    
            vonMises.append( np.sqrt( pow(es[0],2) - es[0]*es[1] + pow(es[1],2) + 3*es[2] ) ) 
            
            ## es: [sigx sigy tauxy]
           # sigmaij = np.array([[es(i,1),es(i,3),0],[es(i,3),es(i,2),0],[0,0,0]])
            sigmaij = np.array([[es[0],es[2],0],
                        [es[2],es[1],0],
                        [0,0,0]])
            [v,w] = np.linalg.eig(sigmaij)
            stresses1.append(v[0]*w[0])
            stresses2.append(v[1]*w[1])
            
        # --- Överför modell variabler till lokala referenser
       
        
        self.outputData.vonMises = vonMises
        self.outputData.edof = edof
        self.outputData.coords = coords
        self.outputData.stresses1 = stresses1
        self.outputData.stresses2 = stresses2
        self.outputData.geometry = geometry
        self.outputData.asolve = asolve
        self.outputData.r = r
        self.outputData.ed = ed
        self.outputData.qs = qs_array
        self.outputData.qt = qt_array
        self.outputData.dofsPerNode = dofsPerNode
        self.outputData.elType = elType
        self.outputData.calcDone = True
    
        
    def executeParamStudy(self):
        """Kör parameter studie"""

        # -- Lagra tidigare värden på d
    
        old_a = self.inputData.a
        old_h = self.inputData.h
        
        i = 1
    
        if self.inputData.paramD:
    
            # --- Skapa värden att simulera över
    
            dRange = np.linspace(self.inputData.dStart, self.inputData.dEnd,
                self.inputData.paramSteps)
    
            # --- Starta parameterstudien
    
            for d in dRange:
                print("Executing for a = %g..." % d)
                filename = "Stresstudy_0"
                # --- Sätt önskad parameter i InputData-instansen
                self.inputData.a = float(d)
                # --- Kör beräkningen 
                solver = Solver(self.inputData,self.outputData)
                solver.execute()
                # --- Exportera vtk-fil
                
                filename += str(i)
                self.exportVtk(filename)
                i = i+1
    
        elif self.inputData.paramT:
            tRange = np.linspace(self.inputData.tStart, self.inputData.tEnd,
                self.inputData.paramSteps)
    
            # --- Starta parameterstudien
    
            for t in tRange:
                print("Executing for h = %g..." % t)
                filename = "Stresstudy_0"
                # --- Sätt önskad parameter i InputData-instansen
                self.inputData.h = float(t)
                # --- Kör beräkningen 
                solver = Solver(self.inputData,self.outputData)
                solver.execute()
                # --- Exportera vtk-fil
                
                filename += str(i)
                self.exportVtk(filename)
                i = i+1
    
        # --- Återställ ursprungsvärden
    
        self.inputData.a = old_a
        self.inputData.h = old_h
    
    def exportVtk(self, filename):
        """Export results to VTK"""
    
        print("Exporting results to %s." % filename)
    
        # --- Skapa punkter och polygon definitioner från vårt nät
    
        points = self.outputData.coords.tolist()
    
        # --- Tänk på att topologin i VTK är 0-baserad varför vi måste minskar **edof** med 1.
    
        #polygons = (self.outputData.edof-1).tolist()
    
        # --- För spänningsproblemet användas, se också nästa stycke:
    
        polygons = (self.outputData.topo-1).tolist()
    
        # --- Resultat från beräkningen skapas i separata objekt. Punkter i vtk.PointData och
        # --- elementdata i vtk.CellData. Nedan anger vi både vektor data och skalärvärden för elementen.
        # --- Tänk på att vektorerna måste ha 3 komponenter, så lägg till detta i beräkningsdelen.
    
        #pointData = vtk.PointData(vtk.Scalars(self.outputData.a.tolist(), name="Nodal Properties"))
        #ellData = vtk.CellData(vtk.Scalars(self.outputData.vonMises, name="Von Mises"), vtk.Vectors(self.outputData.flow, "flow"))
    
        # --- För spänningsproblemet blir det istället (ingen pointData)
    
    
        #HÄR BLIR DET KNAS#
        cellData = vtk.CellData(vtk.Scalars(self.outputData.vonMises, name="mises"), vtk.Vectors(self.outputData.stresses1, "principal stress 1"), vtk.Vectors(self.outputData.stresses2, "principal stress 2"))
    
        # --- Skapa strukturen för elementnätet.
    
        structure = vtk.PolyData(points = points, polygons = polygons)
    
        # --- Lagra allting i en vtk.VtkData instans
    
       # vtkData = vtk.VtkData(structure, pointData, cellData)
    
        # --- För spänningsfallet
    
        vtkData = vtk.VtkData(structure, cellData)
    
        # --- Spara allt till filen
    
        vtkData.tofile(filename, "ascii")
            
        
class Report(object):
    """Klass för presentation av indata och utdata i rapportform."""
    def __init__(self, inputData, outputData):
        self.inputData = inputData
        self.outputData = outputData
        self.report = ""

    def clear(self):
        self.report = ""

    def addText(self, text=""):
        self.report+=str(text)+"\n"

    def __str__(self):
        self.clear()
     
        self.addText()
        self.addText("-------------- Results ----------------------------------")
        self.addText()
        self.addText("Element Displacements:")
        self.addText(ArrayTable(self.outputData.ed, col_width=18, headers=["Element Displacements"]))
        self.addText()
        self.addText("Reaction Forces:")
        self.addText(ArrayTable(self.outputData.r, col_width=20, headers=["Reaction Forces"]))
        self.addText()
        self.addText("Global Displacement:")
        self.addText(ArrayTable(np.asarray(self.outputData.asolve), col_width=20, headers=["Global Displacement"]))
        self.addText()
        self.addText("Strains:")
        self.addText(ArrayTable(np.asarray(self.outputData.qs), col_width=20, headers=["Strains"]))
        self.addText()
        self.addText(ArrayTable(np.asarray(self.outputData.qt), col_width=20, headers=["Strains"]))
        return self.report
    
class ArrayTable(object):
    def __init__(self, arr, col_width=20, decimals=4, headers=[]):
        self.arr = arr
        self.col_width = col_width
        self.decimals = decimals
        self.headers = headers
        self.int_value_format = "|{:>"+str(self.col_width)+"d}"
        self.float_value_format = "|{:>"+str(self.col_width)+"."+str(self.decimals)+"e}"
        self.str_format = "|{:>"+str(self.col_width)+"}"

    def __str__(self):

        rows = self.arr.shape[0]
        cols = self.arr.shape[1]

        ostr = ""

    

        for c in range(cols):
            ostr += "+" + "-" * (self.col_width)
        ostr+="+\n"

        if len(self.headers)>0:
            for c in range(cols):
                if c<len(self.headers):
                    ostr += self.str_format.format(self.headers[c])
                else:
                    ostr += self.str_format.format("")
            ostr+="|\n"

            for c in range(cols):
                ostr += "+" + "-" * (self.col_width)
            ostr+="+\n"

        for r in range(rows):
            for c in range(cols):
                if self.arr.dtype == 'float64':
                    ostr += self.float_value_format.format(self.arr[r, c])
                else:
                    ostr += self.int_value_format.format(self.arr[r, c])

            ostr+="|\n"

        for c in range(cols):
            ostr += "+" + "-" * (self.col_width)
        ostr+="+\n"
       
        return ostr
    
class Visualisation(object):
    def __init__(self, inputData, outputData,calcDone):
        self.inputData = inputData
        self.outputData = outputData
        self.calcDone = calcDone
       # self.calcDone = True


        # --- Variabler som lagrar referenser till öppnade figurer

        self.geomFig = None
        self.meshFig = None
        self.elValueFig = None
        self.nodeValueFig = None
        
   

        
        
    def showGeometry(self):
        
        if self.calcDone == True:
                
            """Visa geometri visualisering"""
    
            geometry = self.outputData.geometry
    
            self.geomFig = cfv.figure(self.geomFig)
            cfv.clf()            
            cfv.drawGeometry(geometry, title="Geometry")
        
        
    def showMesh(self):
         if self.calcDone == True:
             self.meshFig = cfv.figure(self.meshFig) 
             cfv.drawMesh(self.outputData.coords,self.outputData.edof, self.outputData.dofsPerNode,self.outputData.elType, 
             filled=True, title="Mesh") #Draws the mesh.
        
    def showNodalValues(self):
         if self.calcDone == True: 
             self.nodeValueFig = cfv.figure(self.nodeValueFig) 
             cfv.drawDisplacements(self.outputData.asolve, self.outputData.coords, self.outputData.edof, self.outputData.dofsPerNode, self.outputData.elType, 
                          doDrawUndisplacedMesh=False, title="Displacements", 
                          magnfac=25.0)
    def showElementValues(self):
         if self.calcDone ==True:
             self.elValueFig = cfv.figure(self.elValueFig) 
             cfv.drawElementValues(self.outputData.vonMises, self.outputData.coords, self.outputData.edof, self.outputData.dofsPerNode, self.outputData.elType, self.outputData.asolve, 
                          doDrawMesh=True, doDrawUndisplacedMesh=False, 
                          title="Effective Stress", magnfac=25.0)
           
                          
             cfv.colorBar().SetLabel("Effective stress")
    
             cfu.info("Done drawing...")
           # cfv.show()
    def closeAll(self):
        
  
        
        
        self.geomFig = None
        self.meshFig = None
        self.elValueFig = None
        self.nodeValueFig = None
        
        
        cfv.clf()
        cfv.close_all()
        
        
    def wait(self):
        """Denna metod ser till att fönstren hålls uppdaterade och kommer att returnera
        När sista fönstret stängs"""

        cfv.showAndWait()