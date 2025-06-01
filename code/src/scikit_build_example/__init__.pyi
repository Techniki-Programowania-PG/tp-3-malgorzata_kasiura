"""
Pybind11 example plugin
-----------------------

.. currentmodule:: scikit_build_example

.. autosummary::
    :toctree: _generate

    plot_line
    DFT
    IDFT
    splot1d
    sinus
    cosinus
    rectangle
    sawtooth
    filtr1d
    filtr2d
    high_pass_filter
"""

def plot_line(x_axis:  list[float], signal: list[float]) -> None:
    """
    funkcja rysuje wykres  
    """

def DFT(signal: list[float]) -> list[complex]:
    """
    funkcja licz¹ca transformate Fouriera z wektora liczb rzeczywistych
    """

def IDFT(signal: list[complex]) -> list[float]:
    """
    odwrotnoœæ funkcji DFT(), funkcja licz¹ca odwrotn¹ transformate Fouriera do wektora liczb rzeczywistych
    """

def splot1d(x: list[float], h: list[float]) -> list[float]:
    """
    funkcja licz¹ca splot dwóch sygnalow
    """

def sinus(frequency: float, t_start: int, t_end: int, num_samples: int) -> list[float]:
    """
    funkcja zwracaj¹ca sygna³ sinusoidalny   
    """


def cosinus(frequency: float, t_start: int, t_end: int, num_samples: int) -> list[float]:
    """
    funkcja zwracaj¹ca sygna³ cosinusoidalny   
    """

def rectangle(t_start: int, t_end: int, num_samples: int) -> list[float]:
    """
    funkcja zwracaj¹ca sygna³ prostok¹tny   
    """

def sawtooth(frequency: float, t_start: int, t_end: int,num_samples: int) -> list[float]:
    """
    funkcja zwracaj¹ca sygna³ pi³okszta³tny   
    """

def filtr1d(x: list[float],kernel_size:float)->list[float]:
    """
    funkcja zwracaj¹ca przefiltrowany sygna³ jednowymiarowy
    """

def splot2d(x: list[list[float]],h:list[list[float]])->list[list[float]]:
    """
    funkcja licz¹ca splot dwóch macierzy
    """

def filtr2d(x: list[list[float]],kernel_size:float)->list[list[float]]:
    """
    funkcja zwracaj¹ca przefiltrowany sygna³ dwuwymiarowy
    """

def high_pass_filter(signal:  list[float],f0:float,num_samples: int,t_end:float)-> list[float]:
    """
    usuwanie niskich czêstotliwoœci z sygna³u za pomoc¹ DFT
    """

def read_image(nazwa: str)-> int:
    """
    "zamiana obrazu w pliku na macierz wartoœci liczbowych"
    """

def show_image(nazwa: list[list[float]])-> None:
    """
    wyœwietlanie obrazu o zadanej nazwie
    """