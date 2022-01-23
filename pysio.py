# Toolbox Version 0.5.0
import scipy.io
import sys
import getopt
import features as feature
import tkinter as tk
from tkinter import *
import numpy as np
from matplotlib.widgets import SpanSelector
import matplotlib.pyplot as plt
from scipy import signal as sigg
from scipy.signal import find_peaks
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
from tkinter import simpledialog, filedialog
import os

# Initializng parameters to get from the user
filename = ""
dict_name = ""
d_type = ""
fs = 0
x_max = 1
x_min = 0

# Argument imputs
argv = sys.argv[1:]
try:
    opts, args = getopt.getopt(argv, "t:d:n:f:")  
except:
    # More detailed errors will be added on future versions
    print("Error in given arguments")
for opt, arg in opts:
    if opt in ['-t']:
        d_type = arg
        if d_type == "mat":
            print("You chose a .mat file to work on")
        elif d_type == "txt":
            print("You chose a .txt file to work on")
        else:
            print("Wrong data type please enter a valid data type")
    elif opt in ['-d']:
        dict_name = arg
    elif opt in ['-n']:
        filename = arg
    elif opt in ['-f']:
        fs = int(arg)

if filename == "":
    if d_type == "mat":
        if dict_name == "":
            dict_name = input("Please enter the dictionary name of your data: ")
        filetypes = [('all files', '.mat')]
        file_input_window = tk.Tk()
        file_input_window.eval('tk::PlaceWindow . center')
        filename = filedialog.askopenfilename(parent=file_input_window,
                                            initialdir=os.getcwd(),
                                            title="Please select a file:",
                                            filetypes=filetypes)
        file_input_window.destroy()
    else:
        filetypes = [('all files', '.txt')]
        file_input_window = tk.Tk()
        file_input_window.eval('tk::PlaceWindow . center')
        filename = filedialog.askopenfilename(parent=file_input_window,
                                            initialdir=os.getcwd(),
                                            title="Please select a file:",
                                            filetypes=filetypes)
        file_input_window.destroy()
        data = np.loadtxt(filename)
else:
    if filename[-3:] == "mat":
        if dict_name == "":
            dict_name = input("Please enter the dictionary name of your data: ")
        mat = scipy.io.loadmat(filename)
        data = mat.get(dict_name).reshape(-1)
        print("Choosen file:", filename)
    elif filename[-3:] == "txt":
        data = np.loadtxt(filename)
    else:
        raise Exception("Wrong datatype")

if fs == 0:
    # Ask for fs and close window
    fs_input_window = tk.Tk()
    fs_input_window.eval('tk::PlaceWindow . center')
    fs = simpledialog.askinteger("Input", "fs", parent=fs_input_window)
    print("Fs value:", fs)
    fs_input_window.destroy()

window = tk.Tk()
fft = np.zeros((1, 1))
spec = np.zeros((1, 1))
fig = plt.Figure(figsize=(5,4))
region_y = data
dft_frame = tk.Frame(window)
first_dft = True


def select_file():
    filename = tk.filedialog.askopenfilename()

def onselect(xmin, xmax):
    global region_y, x_max, x_min
    xmin /= len(data)
    xmax /= len(data)
    # Discarding the selection beyond data
    if xmin < 0: xmin = 0
    if xmax > 1: xmax = 1
    x_max = xmax
    x_min = xmin
    print("---")
    print("Selected window min:", xmin*len(data))
    print("Selected window max:", xmax*len(data))
    x = np.arange(0.0, 5.0, 0.01)
    indmin, indmax = np.searchsorted(x, (xmin, xmax))
    indmax = min(len(x)-1, indmax)
    indmin -= 1
    indmax -= 1
    # Discarding the selection beyond data
    if indmin < 0: indmin = 0
    if indmax > 99: indmax = 99
    region_y = data[int(indmin*len(data)/99):int(indmax*len(data)/99)]

def create_window(data):
    global ax2, fig

    window.title("PySio - Feature Extraction Tool")
    window.geometry('1200x700')
    canvas = FigureCanvasTkAgg(fig, window)
    canvas.get_tk_widget().grid(row=0,column=0, sticky="nsew")
    ax = fig.add_subplot(111)
    fig.subplots_adjust(bottom=0.25)
    ax.set_xlabel("sample")
    ax.plot(data)
    ax.set_title("Data")
    span = SpanSelector(ax, onselect, 'horizontal', useblit=True,
                        rectprops=dict(alpha=0.5, facecolor='tab:blue'), span_stays=True,)
    plt.show()

    dft() # Calling dft function for seeing the graph at start
    time_domain_text = tk.Label(window, text="Time Domain", font='Helvetica 10 bold')
    time_domain_text.grid(row=1, column=0, sticky="ew")
    spec_domain_text = tk.Label(window, text="Spectral Domain", font='Helvetica 10 bold')
    spec_domain_text.grid(row=6, column=0, sticky="ew")
    Grid.rowconfigure(window, 0, weight=1)
    Grid.columnconfigure(window, 0, weight=1)
    Grid.columnconfigure(window, 1, weight=1)
    Grid.columnconfigure(window, 2, weight=1)
    btn_dft = tk.Button(window, text="DFT", command=dft)
    btn_dft.grid(row=0, column=1, sticky="ew")
    btn1 = tk.Button(window, text="Energy", command=energy)
    btn2 = tk.Button(window, text="ZCR", command=zcr)
    btn3 = tk.Button(window, text="Entropy", command=entropy)
    btn4 = tk.Button(window, text="Time Domain Peaks", command=time_domain_peaks)
    btn5 = tk.Button(window, text="Spectral Entropy", command=lambda:[dft(),spectral_entropy()])
    btn6 = tk.Button(window, text="Spectral Rolloff", command=lambda:[dft(),spectral_rolloff()])
    btn7 = tk.Button(window, text="Spectral Centroid", command=lambda:[dft(),spectral_centroid()])
    btn8 = tk.Button(window, text="Spectral Spread", command=lambda:[dft(),spectral_spread()])
    btn9 = tk.Button(window, text="Spectrogram", command=spectrogram)
    btn10 = tk.Button(window, text="Spectral Peaks", command=lambda:[dft(),spectral_peaks()])
    btn11 = tk.Button(window, text="Bandpower", command=lambda:[dft(),bandpower(data, fs)])
    btn_list = [btn1,btn2,btn3,btn4,btn5,btn6,btn7,btn8,btn9,btn10,btn11]
    row_number = 2
    for btn in btn_list:
        if row_number == 6: row_number += 1
        btn.grid(row=row_number, column=0, sticky="nsew")
        Grid.rowconfigure(window, row_number, weight=1)
        row_number += 1
    window.mainloop()


def energy(): # time domain
    energy = feature.energy(region_y)
    result_text = tk.Text(window, height=1, width=15)
    result_text.grid(row=2, column=1)
    result_text.insert(tk.END, energy)
    result_text.config(state=DISABLED)
    return energy


def zcr(): # time domain
    zcr = feature.zcr(region_y)
    result_text = tk.Text(window, height=1, width=15)
    result_text.grid(row=3, column=1)
    result_text.insert(tk.END, zcr)
    result_text.config(state=DISABLED)
    return zcr


def entropy(): # time domain
    entropy = feature.entropy(region_y)
    result_text = tk.Text(window, height=1, width=15)
    result_text.grid(row=4, column=1)
    result_text.insert(tk.END, entropy)
    result_text.config(state=DISABLED)
    return entropy


def time_domain_peaks():
    result_text = tk.Text(window, height=1, width=15)
    result_text.grid(row=5, column=1)
    result_text.insert(tk.END, "calculated")
    result_text.config(state=DISABLED)
    top_3 = feature.time_domain_peaks(region_y, fs)
    print('Data max 1:', top_3[1], 'Time max 1:', top_3[0] + x_min*len(data))
    print('Data max 2:', top_3[3], 'Time max 2:', top_3[2] + x_min*len(data))
    print('Data max 3:', top_3[5], 'Time max 3:', top_3[4] + x_min*len(data))


def dft(): 
    global fft, first_dft
    fft, f_req = feature.dft(region_y, fs)
    result_text = tk.Text(dft_frame, height=1, width=15)
    if first_dft:
        result_text.pack()
    result_text.insert(tk.END, "FFT calculated")
    result_text.config(state=DISABLED)
    fig2 = plt.Figure(figsize=(5,4))
    canvas2 = FigureCanvasTkAgg(fig2, window)
    canvas2.get_tk_widget().grid(row=0, column=2, sticky="nsew")
    ax2 = fig2.add_subplot(111)
    ax2.set_title("FFT")
    fig2.subplots_adjust(bottom=0.25)
    x = fs*np.arange(0,int(len(region_y)/2),1) / len(region_y)
    ax2.set_xlabel("frequency")
    if len(fft) > len(x):
        fft = fft[:-1]
    ax2.plot(x,fft)
    first_dft = False
    toolbar = NavigationToolbar2Tk(canvas2, window, pack_toolbar=False)
    toolbar.update()
    toolbar.grid(row=1, column=2)
    return fft

def spectrogram():
    global fft
    fft, f_req = feature.dft(region_y, 4000)
    result_text = tk.Text(window, height=1, width=15)
    result_text.grid(row=11, column=1)
    result_text.insert(tk.END, "calculated")
    result_text.config(state=DISABLED)
    winSize = int(np.round(0.1*fs))
    overlap = 90
    noverlap = int(np.round((overlap/100)*winSize))
    f_spect, t_spect, Sxx = sigg.spectrogram(region_y, fs=fs , nperseg=winSize  , noverlap=noverlap,  scaling='spectrum' )
    plt.close()
    fig3 = plt.figure(figsize=(5,4))
    plt.pcolormesh(t_spect, f_spect, 10*np.log10(np.abs(Sxx)) , cmap ='jet', shading='auto')
    plt.ylabel('Frequency [Hz]')
    plt.xlabel('Time [sec]')
    plt.title('Spectrogram')
    plt.colorbar()
    plt.figure(1)
    canvas3 = FigureCanvasTkAgg(fig3, window)
    canvas3.get_tk_widget().grid(row=2, column=2, sticky="nsew", rowspan=13)
    canvas3.draw()
    return fft
    
def spectral_peaks():
    global fft
    fft, f_req = feature.dft(region_y, fs)
    result_text = tk.Text(window, height=1, width=15)
    result_text.grid(row=12, column=1)
    result_text.insert(tk.END, "calculated")
    result_text.config(state=DISABLED)
    peaks, prop = find_peaks(fft, height=0)
    print('Spectral Max Amplitude', np.max(prop['peak_heights']))
    print('Fundamental Frequency:', peaks[np.argmax(prop['peak_heights'])]/2)

def spectral_entropy():
    spectral_entropy = feature.spectral_entropy(fft)
    result_text = tk.Text(window, height=1, width=15)
    result_text.grid(row=7, column=1)
    result_text.insert(tk.END, spectral_entropy)
    result_text.config(state=DISABLED)
    return spectral_entropy


def spectral_rolloff():
    spectral_rolloff = feature.spectral_rolloff(fft)
    result_text = tk.Text(window, height=1, width=15)
    result_text.grid(row=8, column=1)
    result_text.insert(tk.END, spectral_rolloff)
    result_text.config(state=DISABLED)
    return spectral_rolloff


def spectral_centroid():
    c = feature.spectral_centroid(fft, fs)
    result_text = tk.Text(window, height=1, width=15)
    result_text.grid(row=9, column=1)
    result_text.insert(tk.END, c)
    result_text.config(state=DISABLED)
    return spectral_centroid


def spectral_spread():
    s = feature.spectral_spread(fft, fs)
    result_text = tk.Text(window, height=1, width=15)
    result_text.grid(row=10, column=1)
    result_text.insert(tk.END, s)
    result_text.config(state=DISABLED)
    return spectral_spread


def bandpower(x, fs):
    fmin = 0
    fmax = 0
    fmin_window = tk.Tk()
    fmin_window.eval('tk::PlaceWindow . center')
    fmin = int(simpledialog.askinteger("Input", "fmin", parent=fmin_window))
    fmin_window.destroy()
    fmax_window = tk.Tk()
    fmax_window.eval('tk::PlaceWindow . center')
    fmax = int(simpledialog.askinteger("Input", "fmax", parent=fmax_window))
    fmax_window.destroy()   
    f, Pxx = scipy.signal.periodogram(x, fs=fs)
    ind_min = np.argmax(f>fmin) - 1
    ind_max = np.argmax(f>fmax) - 1
    result_text = tk.Text(window, height=1, width=15)
    result_text.grid(row=13, column=1)
    result_text.insert(tk.END, "calculated")
    result_text.config(state=DISABLED)
    print('Bandpower:', np.trapz(Pxx[ind_min:ind_max], f[ind_min:ind_max]))


if __name__ == '__main__':
    create_window(data)
