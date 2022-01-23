# Toolbox Version 0.5.0
import numpy as np
import math
from scipy.signal import argrelextrema, find_peaks


def energy(data):
    return (1 / len(data)) * np.sum(np.square(data))


def zcr(data):
    data_zeros = np.zeros(len(data))
    data_zeros[1:] = data[:len(data) - 1]
    return (1 / (2 * len(data))) * np.sum(np.abs(np.sign(data) - np.sign(data_zeros)))


def entropy(data, num_short_blocks=10):
    eol = np.sum(np.square(data))
    win_len = len(data)
    sub_win_len = math.floor(win_len / num_short_blocks)

    if win_len != sub_win_len * num_short_blocks:
        data = data[0:sub_win_len * num_short_blocks]
    sub_wins = data.reshape(sub_win_len, num_short_blocks, order='F').copy()
    norm_sub_frame_energies = np.zeros((1, sub_wins.shape[1]))
    for i in range(sub_wins.shape[1]):
        norm_sub_frame_energies[0, i] = np.sum(np.square(sub_wins[:, i])) / (eol + np.spacing(1))
    energy_entropy = 0
    for i in range(norm_sub_frame_energies.shape[1]):
        energy_entropy -= norm_sub_frame_energies[0, i] * math.log(norm_sub_frame_energies[0, i] + np.spacing(1), 2)
    return energy_entropy


def time_domain_peaks(data, fs):
    indexes, _ = find_peaks(data)
    index_data_dict = {}
    for i in indexes:
        index_data_dict[i] = data[i]
    sorted_dict = dict(sorted(index_data_dict.items(), key=lambda item: item[1]))
    result = list()
    for i in range(1,4):
        result.append(list(sorted_dict.items())[-1*i][0])
        result.append(list(sorted_dict.items())[-1*i][1])
    return result

def dft(data, f_s=2000, p=0):
    win_len = len(data)
    fft = np.abs(np.fft.fft(data)) / win_len
    if not p:
        fft = fft[0:math.ceil(win_len)]
        f_req = (f_s / 2) * np.arange(0, np.ceil(win_len / 2) + 1) / np.ceil(win_len / 2)
    else:
        fft = np.fft.fftshift(fft)
        if win_len % 2:
            f_req = np.arange(-(win_len - 1) / 2, (win_len - 1) / 2 + 1)
        else:
            f_req = np.arange(-win_len / 2, win_len / 2)
    fft_1 = np.abs(fft)/win_len
    fft_2 = fft_1[1:(round(win_len / 2) + 1)]
    fft_2 = 2*fft_2
    return fft_2, f_req


def spectral_entropy(data, num_short_blocks=10):
    return entropy(data, num_short_blocks)


def spectral_rolloff(data, c=0.90):
    total_energy = np.sum(np.square(data))
    curr_energy = 0
    count_fft = 0
    fft_len = len(data)
    while curr_energy <= c * total_energy and count_fft <= fft_len:
        curr_energy += data[count_fft] ** 2
        count_fft += 1
    count_fft -= 1
    return (count_fft - 1) / fft_len


def spectral_centroid(data, f_s=2000):
    fft_len = len(data)
    m = np.transpose((f_s / (2 * fft_len)) * np.arange(1, fft_len+1))
    data = data / np.max(data)
    c = np.sum(np.multiply(m, data)) / (np.sum(data) + np.spacing(1))
    s = math.sqrt(np.sum(np.square(m - c) * data) / (np.sum(data) + np.spacing(1))) / (f_s / 2)
    c = c / (f_s / 2)
    return c


def spectral_spread(data, f_s=2000):
    fft_len = len(data)
    m = np.transpose((f_s / (2 * fft_len)) * np.arange(1, fft_len+1))
    data = data / np.max(data)
    c = np.sum(np.multiply(m, data)) / (np.sum(data) + np.spacing(1))
    s = math.sqrt(np.sum(np.square(m - c) * data) / (np.sum(data) + np.spacing(1))) / (f_s / 2)
    c = c / (f_s / 2)
    return s
