import scikit_build_example as sbe
import numpy as np

#wizualizacja sygnalu 1d
x = np.linspace(0, 10 * np.pi, 1000)
y = np.sin(x)
#sbe.plot_line(x, y)

#DFT
f = 10
N = 1000
t = np.linspace(0, 2 * np.pi, N)
x = np.sin(f * t)
X = np.abs(sbe.DFT(x)) ** 2
#sbe.plot_line(t, x)
# zapisane jako sin_ft.png
#sbe.plot_line(np.arange(N), X)
# zapisane jako sin_ft_fft.png
x = np.sin(10 * f * t)
X = np.abs(sbe.DFT(x)) ** 2
#sbe.plot_line(t, x)
# zapisane jako sin_10_ft.png
#sbe.plot_line(np.arange(N), X)
# zapisane jako sin_10_ft_fft.png

#IDFT
f = 10
N = 1000
t = np.linspace(0, 2 * np.pi, N)
x = np.sin(f * t)
x1=np.abs(sbe.DFT(x))**2
x2 = sbe.IDFT(sbe.DFT(x))
#sbe.plot_line(t,x)
# zapisane jako sin_ft.png
#sbe.plot_line(t,x1)
# zapisane jako sin_ft_fft.png
#sbe.plot_line(t,x2)
# zapisane jako sin_10_ft_ifft.png

x = np.sin(10 * f * t)
x1 = np.abs(sbe.DFT(x)) ** 2
x2=sbe.IDFT(sbe.DFT(x))
#sbe.plot_line(t, x)
# zapisane jako sin_10_ft.png
#sbe.plot_line(np.arange(N), x1)
# zapisane jako sin_10_ft_fft.png
#sbe.plot_line(np.arange(N), x2)
# zapisane jako sin_10_ft_ifft.png

#filtr 1d
step_signal = np.zeros(100)
step_signal[50:] = 1
x=np.linspace(0,1,100)
#sbe.plot_line(x,step_signal)
# zapisane jako step_signal.png
noisy_signal = (step_signal + np.random.normal(0, 0.33,len(step_signal)))
#sbe.plot_line(x,noisy_signal)
# zapisane jako noisy_signal.png
filtered_signal_3=sbe.filtr1d(noisy_signal,3)
#sbe.plot_line(x,filtered_signal_3)
# zapisane jako filtered_signal_3.png
filtered_signal_30=sbe.filtr1d(noisy_signal,30)
#sbe.plot_line(x,filtered_signal_30)
# zapisane jako filtered_signal_30.png

#filtr 2d
x=sbe.read_image("kot.jpeg")
#sbe.show_image(sbe.read_image("kot.jpeg"))
# zapisane jako image.png
filtered_image_2=sbe.filtr2d(x,2)
#sbe.show_image(filtered_image_2)
# zapisane jako filtered_image_2.png
filtered_image_10=sbe.filtr2d(x,10)
#sbe.show_image(filtered_image_10)
# zapisane jako filtered_image_10.png

#sinus
y1=sbe.sinus(1/15,15,48,1000)
x1=np.linspace(0,45,len(y1))
#sbe.plot_line(x1,y1)
# zapisane jako sinus.png

y2=sbe.cosinus(1.5,2,4,1000)
x2=np.linspace(0,4,len(y2))
#sbe.plot_line(x2,y2)
# zapisane jako cosinus.png

y3=sbe.rectangle(8,17,1000)
x3=np.linspace(0,17,len(y3))
#sbe.plot_line(x3,y3)
# zapisane jako rectangle.png

y4=sbe.sawtooth(1/6,4,26,1000)
x4=np.linspace(0,26,len(y4))
#sbe.plot_line(x4,y4)
# zapisane jako sawtooth.png

#filtr gornoprzepustowy
s1=sbe.sinus(2,0,50,1000)
s2=sbe.sinus(1/10,0,50,1000)
y=np.zeros(len(s1))
for i in range(len(s1)):
    y[i]=0.5*s1[i]+s2[i]

x=np.linspace(0,50,len(y))
y1=np.abs(sbe.DFT(y))**2
sbe.plot_line(x,y)
# zapisane jako suma_sinus.png
sbe.plot_line(x,y1)
# zapisane jako dft_suma_sinus.png
y2=sbe.high_pass_filter(y,2,1000,50)
# zapisane jako fir_suma_sinus.png
sbe.plot_line(x,y2)
# zapisane jako dft_firsuma_sinus.png

#wykrywanie krawêdzi
