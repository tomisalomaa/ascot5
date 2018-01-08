"""
ASCOT5 graphical interface.
"""
import matplotlib
matplotlib.use("TkAgg")
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2TkAgg
from matplotlib.figure import Figure

import tkinter as tk

import numpy as np

import a5py.ascot5io.ascot5 as ascot5
import a5py.postprocessing.state as stateplot

# The GUI itself does not plot anything - it should be done 
# in these imported packages.
import a5py.gui.guistate as guistate
import a5py.gui.guiplasma as guiplasma
import a5py.gui.guibfield as guibfield

LARGE_FONT= ("Verdana", 12)


class ascot5GUI(tk.Tk):
    """
    The main window. During initialization all main frames are initialized.
    They are:
    
    Menu     - The initial view
    
    Input    - For input data overview
    Plasma   - Plasma input
    Bfield   - Magnetic field input
    
    Inistate - Displays simulation inistate
    Endstate - Displays simulation endstate
    
    Additionally we use a frame dubbed WIP to indicate Work In Progress features.
    """
    
    def __init__(self, a5fn, *args, **kwargs):
        """
        Initialize the main window.
        
        Parameters
        ----------
        
        a5fn : str
            Full path to ASCOT5 HDF5 file this GUI corresponds to.        
        """
        tk.Tk.__init__(self, *args, **kwargs)
        tk.Tk.wm_title(self, "ASCOT5 GUI")

        # Specify window size in pixels
        self.geometry('1200x800+0+0')
    
        # The full path to an ASCOT5 HDF5 file we are dealing with.
        self.a5fn = a5fn

        # The frame which displays all other frames.
        container = tk.Frame(self)
        container.pack(side="top", fill="both", expand = True)
        container.grid_rowconfigure(0, weight=1)
        container.grid_columnconfigure(0, weight=1)

        # Initialize all frames.
        self.frames = {}
        for F in (FrameMenu, FrameWIP, FrameInistate, FrameEndstate, FrameInput,\
                  FramePlasma, FrameBfield):

            frame = F(container, self)

            self.frames[F] = frame

            frame.grid(row=0, column=0, sticky="nsew")

        self.show_frame(FrameMenu)

    def show_frame(self, cont):
        """
        Display the given frame.
        
        Parameters
        ----------
        
        cont : Frame object
            The Frame to be displayed.
        """
        frame = self.frames[cont]
        frame.tkraise()

        
class FrameMenu(tk.Frame):
    """
    The navigation frame.
    
    TODO: Overlay is an ugly mess (functionality is okay)
    """
    def __init__(self, parent, controller):
        tk.Frame.__init__(self,parent)
        label = tk.Label(self, text="Main", font=LARGE_FONT)
        label.pack(pady=10,padx=10)


        # Divide this frame first in to two frames vertically and then divide
        # both horizontally to three frames to get a total of six sub-frames.
        frameN = tk.Frame(self)
        frameN.pack(side = "top", fill = "both", expand=True)

        frameS = tk.Frame(self)
        frameS.pack(side = "top", fill = "both", expand=True)
 
        frameNW = tk.Frame(frameN)
        frameNW.pack(side = "left", fill="both", expand=1)
        frameNM = tk.Frame(frameN)
        frameNM.pack(side = "left", fill="both", expand=1)
        frameNE = tk.Frame(frameN)
        frameNE.pack(side = "left", fill="both", expand=1)

        frameSW = tk.Frame(frameS)
        frameSW.pack(side = "left", fill="both", expand=1)
        frameSM = tk.Frame(frameS)
        frameSM.pack(side = "left", fill="both", expand=1)
        frameSE = tk.Frame(frameS)
        frameSE.pack(side = "left", fill="both", expand=1)

        # The North-West region is for choosing the HDF5 file and ASCOT5 binary.
        label = tk.Label(frameNW, text="General", font=LARGE_FONT)
        label.grid(column=0,row=0)
        
        label = tk.Label(frameNW, text="HDF5 file:")
        label.grid(column=0,row=1)
        entry = tk.Entry(frameNW)
        entry.grid(column=1,row=1)
        label = tk.Label(frameNW, text="ASCOT5 binary:")
        label.grid(column=0,row=2)
        entry = tk.Entry(frameNW)
        entry.grid(column=1,row=2)

        # The North region is for input handling.
        label = tk.Label(frameNM, text="Input configuration", font=LARGE_FONT)
        label.pack(expand=1)
        button = tk.Button(frameNM, text="Options",
                            command=lambda: controller.show_frame(FrameWIP))
        button.pack(side = "top",fill='x')
        button = tk.Button(frameNM, text="Choose inputs",
                            command=lambda: controller.show_frame(FrameWIP))
        button.pack(side = "top",fill='x')
        button = tk.Button(frameNM, text="Examine inputs",
                            command=lambda: controller.show_frame(FrameInput))
        button.pack(side = "top",fill='x')
        button = tk.Button(frameNM, text="Add / remove inputs",
                            command=lambda: controller.show_frame(FrameWIP))
        button.pack(side = "top",fill='x')
        
        # The North-East region is for sanity checks and code testing.
        label = tk.Label(frameNE, text="Other", font=LARGE_FONT)
        label.pack(expand=1)
        button = tk.Button(frameNE, text="Sanity checks",
                            command=lambda: controller.show_frame(FrameWIP))
        button.pack(side = "top",fill='x')
        button = tk.Button(frameNE, text="Test ascot",
                            command=lambda: controller.show_frame(FrameWIP))
        button.pack(side = "top",fill='x')

        # The South-West region is for pre-defined runs.
        label = tk.Label(frameSW, text="Pre-defined runs", font=LARGE_FONT)
        label.pack(expand=1)
        button = tk.Button(frameSW, text="Slowing-down simulation",
                            command=lambda: controller.show_frame(FrameWIP))
        button.pack(side = "top",fill='x')
        button = tk.Button(frameSW, text="Transport evaluation",
                            command=lambda: controller.show_frame(FrameWIP))
        button.pack(side = "top",fill='x')
        button = tk.Button(frameSW, text="Poincare plots",
                            command=lambda: controller.show_frame(FrameWIP))
        button.pack(side = "top",fill='x')

        # The South region is for showing the results.
        label = tk.Label(frameSM, text="Results", font=LARGE_FONT)
        label.pack(expand=1)
        button = tk.Button(frameSM, text="Inistate",
                            command=lambda: controller.show_frame(FrameInistate))
        button.pack(side = "top",fill='x')
        button = tk.Button(frameSM, text="Endstate",
                            command=lambda: controller.show_frame(FrameEndstate))
        button.pack(side = "top",fill='x')
        button = tk.Button(frameSM, text="Distributions",
                            command=lambda: controller.show_frame(FrameWIP))
        button.pack(side = "top",fill='x')
        button = tk.Button(frameSM, text="Orbits",
                            command=lambda: controller.show_frame(FrameWIP))
        button.pack(side = "top",fill='x')

        # The South-East region is for postprocessing.
        label = tk.Label(frameSE, text="Postprocessing", font=LARGE_FONT)
        label.pack(expand=1)
        button = tk.Button(frameSE, text="Add distributions",
                            command=lambda: controller.show_frame(FrameWIP))
        button.pack(side = "top",fill='x')
        button = tk.Button(frameSE, text="Add state fields",
                            command=lambda: controller.show_frame(FrameWIP))
        button.pack(side = "top",fill='x')


class FrameWIP(tk.Frame):
    """
    Temporary frame indicating work in progress
    """

    def __init__(self, parent, controller):
        tk.Frame.__init__(self, parent)
        label = tk.Label(self, text="This feature has not been implemented yet", font=LARGE_FONT)
        label.pack(pady=10,padx=10)

        button = tk.Button(self, text="Back to Main",
                            command=lambda: controller.show_frame(FrameMenu))
        button.pack()


class FrameState(tk.Frame):
    """
    Frame for displaying marker state.
    """

    def __init__(self, parent, controller, mrkstate):
        tk.Frame.__init__(self, parent)
        self.mrkstate = mrkstate
        self.a5fn = controller.a5fn
        self.plotfig = Figure(figsize=(5,5), dpi=100)

        self.showparticles = False
        self.showlogx      = False
        self.showlogy      = False

        # Sometimes one likes to see marker inistate versus marker endCond, this flag enables that
        self.notendstate = self.mrkstate != "endstate"
        self.showendcond = False
        
        # Build labels and buttons
        label = tk.Label(self, text=self.mrkstate, font=LARGE_FONT)
        label.pack(pady=10,padx=10)

        button = tk.Button(self, text="Back to Main",
                           command=lambda: controller.show_frame(FrameMenu))
        button.pack()

        def toggleparticles():
            self.showparticles=(not self.showparticles)
            self.tkraise()
        button = tk.Button(self, text="Show particles",
                           command=lambda: toggleparticles() )
        button.pack()

        def togglelogx():
            self.showlogx=(not self.showlogx)
            self.tkraise()
        button = tk.Button(self, text="Logarithmic x-axis",
                           command=lambda: togglelogx() )
        button.pack()

        def togglelogy():
            self.showlogy=(not self.showlogy)
            self.tkraise()
        button = tk.Button(self, text="Logarithmic y-axis",
                           command=lambda: togglelogy() )
        button.pack()

        if self.notendstate:
            def toggleendcond():
                self.showendcond=(not self.showendcond)
                self.tkraise()
            button = tk.Button(self, text="Show endcond",
                               command=lambda: toggleendcond() )
            button.pack()

        # Build a canvas and add a tool bar to it
        canvas = FigureCanvasTkAgg(self.plotfig, self)
        canvas.get_tk_widget().pack(side=tk.BOTTOM, fill=tk.BOTH, expand=True)

        toolbar = NavigationToolbar2TkAgg(canvas, self)
        toolbar.update()
        canvas._tkcanvas.pack(side=tk.TOP, fill=tk.BOTH, expand=True)

        self.canvas = canvas

    def tkraise(self):
        super(FrameState, self).tkraise()

        mrkendcond = []

        if self.showendcond:
            s = ascot5.read_hdf5(self.a5fn, "states")["states"]["endstate"]
            mrkendcond = s["endCond"]

        s = ascot5.read_hdf5(self.a5fn, "states")["states"][self.mrkstate]
        stateplot.gather(s)
        
        if not self.showendcond:
            mrkendcond = s["endCond"]
            
        eid = np.unique(mrkendcond)
        endcond = []; Ekin = []; pitch = []; time = []; rho = []; cputime = []; phi = []; charge = []; weights = []; R = []; z = []
        endcondlabel = []
        for i in eid:
            ix = np.argwhere(np.isin(mrkendcond, i))
            endcond.append(np.take(mrkendcond,ix))
            Ekin.append(np.take(s["Ekin"],ix))
            pitch.append(np.take(s["pitch"],ix))
            time.append(np.take(s["time"],ix))
            rho.append(np.take(s["rho"],ix))
            cputime.append(np.take(s["cputime"],ix))
            phi.append(np.take(s["phi"],ix))
            charge.append(np.take(s["charge"],ix))
            weights.append(np.take(s["weight"],ix))
            R.append(np.take(s["R"],ix))
            z.append(np.take(s["z"],ix))

            endcondlabel.append(stateplot.interpEndcond(i))

        if not self.showparticles:
            weights = None
            
        if self.showlogx:
            for i in range(len(Ekin)):
                Ekin[i] = np.log10(Ekin[i])
            for i in range(len(time)):
                time[i] = np.log10(time[i]+1e-12)
            for i in range(len(cputime)):
                cputime[i] = np.log10(cputime[i]+1e-12)

                
                
        f = self.plotfig
        f.clf()
        a = f.add_subplot(3,4,1)
        guistate.plotendcond(endcond, endcondlabel, ax=a, logy=self.showlogy, weights=weights)
        
        a = f.add_subplot(3,4,2)
        guistate.plottime(time, ax=a, logy=self.showlogy, weights=weights)
        
        #a = f.add_subplot(3,4,3)

        a = f.add_subplot(3,4,5)
        guistate.plotEkin(Ekin, ax=a, logy=self.showlogy, weights=weights)
        
        a = f.add_subplot(3,4,6)
        guistate.plotpitch(pitch, ax=a, logy=self.showlogy, weights=weights)
        
        #a = f.add_subplot(3,4,7)

        a = f.add_subplot(3,4,9)
        guistate.plotrho(rho, ax=a, logy=self.showlogy, weights=weights)
        
        a = f.add_subplot(3,4,10)
        guistate.plotphi(phi, ax=a, logy=self.showlogy, weights=weights)
        
        a = f.add_subplot(3,4,11)
        guistate.plotcharge(charge, ax=a, logy=self.showlogy, weights=weights)
        
        a = f.add_subplot(3,4,12)
        guistate.plotcputime(cputime, ax=a, logy=self.showlogy, weights=weights)

        a = f.add_subplot(3,4,(4,8))
        guistate.plotRz(R,z,ax=a)
        
        f.tight_layout()
        
        self.canvas.show()


class FrameInistate(FrameState):

    def __init__(self, parent, controller):
        FrameState.__init__(self, parent, controller, "inistate")


class FrameEndstate(FrameState):

    def __init__(self, parent, controller):
        FrameState.__init__(self, parent, controller, "endstate")


class FrameInput(tk.Frame):
    """
    Frame for input setup and analysis.
    """

    def __init__(self, parent, controller):
        tk.Frame.__init__(self, parent)
        label = tk.Label(self, text="Inputs", font=LARGE_FONT)
        label.pack(pady=10,padx=10)

        button = tk.Button(self, text="Back to Main",
                            command=lambda: controller.show_frame(FrameMenu))
        button.pack()
        
        button = tk.Button(self, text="Magnetic field",
                            command=lambda: controller.show_frame(FrameBfield))
        button.pack()
        
        button = tk.Button(self, text="Electric field",
                            command=lambda: controller.show_frame(FrameWIP))
        button.pack()
        
        button = tk.Button(self, text="Plasma",
                            command=lambda: controller.show_frame(FramePlasma))
        button.pack()
        
        button = tk.Button(self, text="Wall",
                            command=lambda: controller.show_frame(FrameWIP))
        button.pack()
      
      
class FramePlasma(tk.Frame):
    """
    Frame for plasma input setup and analysis.
    """

    def __init__(self, parent, controller):
        tk.Frame.__init__(self, parent)
        label = tk.Label(self, text="Plasma", font=LARGE_FONT)
        label.pack(pady=10,padx=10)

        button = tk.Button(self, text="Back to Main",
                            command=lambda: controller.show_frame(FrameMenu))
        button.pack()
        
        button = tk.Button(self, text="Back to Inputs",
                            command=lambda: controller.show_frame(FrameInput))
        button.pack()
        
        self.a5fn = controller.a5fn
        self.plotfig = Figure(figsize=(5,5), dpi=100)
        
        # Build a canvas and add a tool bar to it
        canvas = FigureCanvasTkAgg(self.plotfig, self)
        canvas.get_tk_widget().pack(side=tk.BOTTOM, fill=tk.BOTH, expand=True)

        toolbar = NavigationToolbar2TkAgg(canvas, self)
        toolbar.update()
        canvas._tkcanvas.pack(side=tk.TOP, fill=tk.BOTH, expand=True)

        self.canvas = canvas

    def tkraise(self):
        super(FramePlasma, self).tkraise()
        pls = ascot5.read_hdf5(self.a5fn, "plasma")["plasma"]
        
        labels_temperature = ["electron"]
        
        labels_density = []
        for i in range(pls["Nion"]):
            labels_density.append("q = " + str(pls["Z_num"][i]) + ", m = " + str(pls["A_mass"][i]))
        
        f = self.plotfig
        f.clf()
        a = f.add_subplot(2,1,1)
        guiplasma.plottemperature(pls["rho"], pls["etemp"], labels_temperature, ax=a)
        
        a = f.add_subplot(2,1,2)
        guiplasma.plotdensity(pls["rho"], pls["idens"], labels_density, ax=a)
        
        f.tight_layout()
        
        self.canvas.show()
        
class FrameBfield(tk.Frame):
    """
    Frame for plasma input setup and analysis.
    """

    def __init__(self, parent, controller):
        tk.Frame.__init__(self, parent)
        label = tk.Label(self, text="Magnetic field", font=LARGE_FONT)
        label.pack(pady=10,padx=10)

        button = tk.Button(self, text="Back to Main",
                            command=lambda: controller.show_frame(FrameMenu))
        button.pack()
        
        button = tk.Button(self, text="Back to Inputs",
                            command=lambda: controller.show_frame(FrameInput))
        button.pack()
        
        self.a5fn = controller.a5fn
        self.plotfig = Figure(figsize=(5,5), dpi=100)
        
        # Build a canvas and add a tool bar to it
        canvas = FigureCanvasTkAgg(self.plotfig, self)
        canvas.get_tk_widget().pack(side=tk.BOTTOM, fill=tk.BOTH, expand=True)

        toolbar = NavigationToolbar2TkAgg(canvas, self)
        toolbar.update()
        canvas._tkcanvas.pack(side=tk.TOP, fill=tk.BOTH, expand=True)

        self.canvas = canvas

    def tkraise(self):
        super(FrameBfield, self).tkraise()
        bfield = ascot5.read_hdf5(self.a5fn, "bfield")["bfield"]["B_2D"]
        
        R = np.linspace(bfield["Rmin"], bfield["Rmax"], bfield["nR"])
        z = np.linspace(bfield["zmin"], bfield["zmax"], bfield["nz"])
        
        rho = np.sqrt( np.abs( (bfield["psi"] - bfield["psiaxis"]) / (bfield["psiaxis"] -  bfield["psisepx"]) ) )
        
        f = self.plotfig
        f.clf()
        
        a = f.add_subplot(2,3,1)
        guibfield.plotsurfandcontour(R, z, bfield["psi"], 'Poloidal flux', ax=a, fig=f)
        
        a = f.add_subplot(2,3,2)
        guibfield.plotsurfandcontour(R, z, rho, 'Normalized poloidal flux', ax=a, fig=f)
        
        a = f.add_subplot(2,3,4)
        guibfield.plotsurf(R, z, bfield["B_R"], 'B_R', ax=a, fig=f)
        
        a = f.add_subplot(2,3,5)
        guibfield.plotsurf(R, z, bfield["B_phi"], 'B_phi', ax=a, fig=f)
        
        a = f.add_subplot(2,3,6)
        guibfield.plotsurf(R, z, bfield["B_z"], 'B_z', ax=a, fig=f)
        
        f.tight_layout()
        
        self.canvas.show()
        

def run(a5fn):
    # Launch GUI
    app = ascot5GUI(a5fn)
    app.mainloop()
