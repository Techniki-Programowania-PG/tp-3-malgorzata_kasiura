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
    funkcja licz�ca transformate Fouriera z wektora liczb rzeczywistych
    """

def IDFT(signal: list[complex]) -> list[float]:
    """
    odwrotno�� funkcji DFT(), funkcja licz�ca odwrotn� transformate Fouriera do wektora liczb rzeczywistych
    """

def splot1d(x: list[float], h: list[float]) -> list[float]:
    """
    funkcja licz�ca splot dw�ch sygnalow
    """

def sinus(frequency: float, t_start: int, t_end: int, num_samples: int) -> list[float]:
    """
    funkcja zwracaj�ca sygna� sinusoidalny   
    """


def cosinus(frequency: float, t_start: int, t_end: int, num_samples: int) -> list[float]:
    """
    funkcja zwracaj�ca sygna� cosinusoidalny   
    """

def rectangle(t_start: int, t_end: int, num_samples: int) -> list[float]:
    """
    funkcja zwracaj�ca sygna� prostok�tny   
    """

def sawtooth(frequency: float, t_start: int, t_end: int,num_samples: int) -> list[float]:
    """
    funkcja zwracaj�ca sygna� pi�okszta�tny   
    """

def filtr1d(x: list[float],kernel_size:float)->list[float]:
    """
    funkcja zwracaj�ca przefiltrowany sygna� jednowymiarowy
    """

def splot2d(x: list[list[float]],h:list[list[float]])->list[list[float]]:
    """
    funkcja licz�ca splot dw�ch macierzy
    """

def filtr2d(x: list[list[float]],kernel_size:float)->list[list[float]]:
    """
    funkcja zwracaj�ca przefiltrowany sygna� dwuwymiarowy
    """

def high_pass_filter(signal:  list[float],f0:float,num_samples: int,t_end:float)-> list[float]:
    """
    usuwanie niskich cz�stotliwo�ci z sygna�u za pomoc� DFT
    """

def read_image(nazwa: str)-> int:
    """
    "zamiana obrazu w pliku na macierz warto�ci liczbowych"
    """

def show_image(nazwa: list[list[float]])-> None:
    """
    wy�wietlanie obrazu o zadanej nazwie
    """