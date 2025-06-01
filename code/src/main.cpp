#include <pybind11/pybind11.h>
#include <matplot\matplot.h>
#include <vector>
#include <pybind11/stl.h>
#include <complex>
#include <pybind11/complex.h>
#include <pybind11/functional.h>
#include <string>

#define STRINGIFY(x) #x
#define MACRO_STRINGIFY(x) STRINGIFY(x)

namespace py = pybind11;


void plot_line(std::vector<double> x_axis, std::vector<double> signal) {

    using namespace matplot;

    plot(x_axis, signal);
    show();
}

std::vector<std::complex<double>> DFT(std::vector<double> signal) {

    using namespace matplot;

    std::vector<std::complex<double>> y;
    std::vector<std::complex<double>> y2;

    y.reserve(signal.size());

    for (int i = 0; i < signal.size(); ++i) {
        std::complex<double> iNum(signal[i], 0);
        y.push_back(iNum);
    }
   
    for (int k = 0; k < signal.size(); ++k) {
        double re = 0.0, im = 0.0;
        for (int n = 0; n < signal.size(); ++n) {
            re += signal[n] * cos(2*pi*k*n/ signal.size());
        }
        for (int n = 0; n < signal.size(); ++n) {
            im -= signal[n] * sin(2 * pi * k * n / signal.size());
        }
        std::complex<double> iNum(re, im);
        y[k]=iNum;
    }
    for (int l = 0; l < signal.size()/2; ++l) {
        y2.push_back(y[l]);
    }
    return y2;
}

std::vector<double> IDFT(std::vector<std::complex<double>> signal) {
    using namespace matplot;
    using namespace std::complex_literals;

    std::vector<double> y;
    y.reserve(signal.size());

    for (int i = 0; i < signal.size(); ++i) {
        y.push_back(0);
    }

    float N = y.size();
    for (int n = 0; n < signal.size(); ++n) {
        double re = 0.0, im = 0.0;
        for (int k = 0; k < signal.size(); ++k) {
            re += signal[k].real() * cos(2 * pi * k * n / signal.size()) - signal[k].imag() * sin(2 * pi * k * n / signal.size());
        }
        y[n] = re/N;
    }
    return y;
}

std::vector<double> splot1d(std::vector<double> x, std::vector<double> h) {
    
    std::vector<double> y;
    for (int i = 0; i < x.size()+h.size()-1; ++i) {
        y.push_back(0);
    }

    for (int n = 0; n < x.size()+h.size()-1; ++n) {
        for (int k = 0; k < x.size(); ++k) {
            if ( n - k < h.size()) {
                y[n] += x[k] * h[n-k];
            }
        }
    }

    return y;
}


std::vector<double> sinus(float frequency, int t_start, int t_end, int num_samples) {
    using namespace matplot;
    double T = 1 / frequency;

    std::vector<double> x = linspace(0, t_end, num_samples);

    std::vector<double> y;
    for (int i = 0; i < x.size(); ++i) {
        y.push_back(0);
    }

    for (int i = 0; i < x.size(); ++i) {
        if (x[i] >= t_start) {
            double temp = x[i] * (2 * pi) / (T);
            y[i] = sin(temp);
        }
    }

    return y;
}

std::vector<double> cosinus(float frequency, int t_start, int t_end, int num_samples) {
    using namespace matplot;
    double T = 1 / frequency;

    std::vector<double> x = linspace(0, t_end, num_samples);

    std::vector<double> y;
    for (int i = 0; i < x.size(); ++i) {
        y.push_back(0);
    }

    for (int i = 0; i < x.size(); ++i) {
        if (x[i] >= t_start) {
            double temp = x[i] * (2 * pi) / (T);
            y[i] = cos(temp);
        }
    }
 
    return y;
}

std::vector<double> rectangle( int t_start, int t_end, int num_samples) {
    using namespace matplot;

    std::vector<double> x = linspace(0, t_end, num_samples);

    std::vector<double> y;
    for (int i = 0; i < x.size(); ++i) {
        y.push_back(0);
    }

    for (int i = 0; i < x.size(); ++i) {
        if (x[i] >= t_start) {
            y[i] = 1;
        }
    }

    return y;
}

std::vector<double> sawtooth(float frequency,int t_start, int t_end, int num_samples) {
    using namespace matplot;
    double T = 1 / frequency;

    std::vector<double> x = linspace(0, t_end, num_samples);

    std::vector<double> y;
    for (int i = 0; i < (t_start*num_samples)/(t_end); ++i) {
        y.push_back(0);
    }

    for (int n = 0; n < x.size()- (t_start * num_samples) / (t_end); n+=(T*num_samples)/(t_end)) {
        for (int i = 0; i < (T * num_samples) / (t_end)-1; i++) {
            y.push_back(x[i]*frequency);
        }
    }

;
    return y;
}

std::vector<double> filtr1d(std::vector<double> x,float a) {
    std::vector<double> y;
    std::vector<double> z;
    for (int i = 0; i < a ; ++i) {
        y.push_back(1/a);
    }

    z = splot1d(x, y);

    return z;

}

std::vector<std::vector<double>> splot2d(std::vector<std::vector<double>> x, std::vector<std::vector<double>> h) {
    std::vector<std::vector<double>> y;


    for (int i = 0; i < x.size() + h.size() - 1; ++i) {
        std::vector<double> temp;
        for (int i = 0; i < x.size() + h.size() - 1; ++i) {
            temp.push_back(0);
        }
        y.push_back(temp);
    }

    for (int n = 0; n < x.size() + h.size() - 1; ++n) {
        for (int k = 0; k < x.size(); ++k) {
            if (n - k < h.size()) {
                for (int i = 0; i < y[n].size(); i++) {
                    y[n][i] += splot1d(h[n - k], x[k])[i];
                }
            }
        }
    }

    return y;
}

std::vector<std::vector<double>> filtr2d(std::vector<std::vector<double>> x, double a){
    std::vector<std::vector<double>> y;
    std::vector<std::vector<double>> z;

    for (int i = 0; i < a; ++i) {
        std::vector<double> temp;
        for (int i = 0; i < a; ++i) {
            temp.push_back(1/a);
        }
        y.push_back(temp);
    }

    z = splot2d(x, y);

    return z;

}



std::vector <double> high_pass_filter(std::vector <double> signal, float f0,int num_samples,float t_end){
    using namespace matplot;
    std::vector<std::complex<double>> z;
    std::complex<double> temp = (0.0,  0.0);    
    std::vector<std::complex<double>> y = DFT(signal);

    for (int i; i < signal.size(); i++) {
        if (    sqrt(pow(y[i].real(),2) + pow(y[i].imag(),2))  <((1/f0 * num_samples) / (t_end))) {
            z.push_back(temp);
        }
        else {
            z.push_back(y[i]);
        }
    }
    std::vector<double> x = IDFT(z);
    return x;
}

std::vector<std::vector<double>> read_image(std::string nazwa) {
    using namespace matplot;
    auto image = imread(nazwa);
    auto x = image[0];
   // imshow(x);
  //  show();
   // return x;

    std::vector<std::vector<double>> y;
    for (int i = 0; i < x.size(); i++) {
        std::vector<double> temp = {};
        for (int j = 0; j < x[i].size(); j++) {
            temp.push_back(x[i][j]);
        }
        y.push_back(temp);
    }
    return y;

}

void show_image(std::vector<std::vector<double>> x) {
    using namespace matplot;


    image(x);

    show();
}


PYBIND11_MODULE(_core, m) {
    m.doc() = "pybind11 example plugin";
    m.def("plot_line", &plot_line,"funkcja rysuje wykres ");
    m.def("DFT", &DFT,"funkcja licz¹ca transformate Fouriera z wektora liczb rzeczywistych");
    m.def("IDFT", &IDFT,"odwrotnoœæ funkcji DFT(), funkcja licz¹ca odwrotn¹ transformate Fouriera do wektora liczb rzeczywistych");
    m.def("splot1d", &splot1d,"funkcja licz¹ca splot dwóch sygnalow");
    m.def("sinus", &sinus,"funkcja zwracaj¹ca sygna³ sinusoidalny");
    m.def("cosinus", &cosinus,"funkcja zwracaj¹ca sygna³ cosinusoidalny");
    m.def("rectangle", &rectangle,"funkcja zwracaj¹ca sygna³ prostok¹tny");
    m.def("sawtooth", &sawtooth,"funkcja zwracaj¹ca sygna³ pi³okszta³tny");
    m.def("filtr1d", &filtr1d,"funkcja zwracaj¹ca przefiltrowany sygna³ jednowymiarowy");
    m.def("splot2d", &splot2d,"funkcja licz¹ca splot dwóch macierzy");
    m.def("filtr2d", &filtr2d,"funkcja zwracaj¹ca przefiltrowany sygna³ dwuwymiarowy");
    m.def("high_pass_filter", &high_pass_filter,"usuwanie niskich czêstotliwoœci z sygna³u za pomoc¹ DFT");
    m.def("read_image", &read_image,"zamiana obrazu w pliku na macierz wartoœci liczbowych");
    m.def("show_image", &show_image,"wyœwietlanie obrazu o zadanej nazwie");


#ifdef VERSION_INFO
    m.attr("__version__") = MACRO_STRINGIFY(VERSION_INFO);
#else
    m.attr("__version__") = "dev";
#endif
}
