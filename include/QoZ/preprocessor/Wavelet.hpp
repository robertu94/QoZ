#ifndef SZ3_WAVELET_HPP
#define SZ3_WAVELET_HPP

#ifdef ENABLE_GSL 
#include "QoZ/preprocessor/PreProcessor.hpp"
#include "QoZ/preprocessor/CDF97.h"
#include <gsl/gsl_wavelet.h>



namespace QoZ {
    template<class T, uint N>

    class Wavelet : public concepts::PreprocessorInterface<T, N> {
    public:
        

        
        void preProcess(T *data, size_t n) {




            
            size_t m = n - 1;
            m |= m >> 1;
            m |= m >> 2;
            m |= m >> 4;
            m |= m >> 8;
            m |= m >> 16;
            m++;

            std::vector<double> dwtdata(m, 0);
            gsl_wavelet *w;
            gsl_wavelet_workspace *work;

            w = gsl_wavelet_alloc(gsl_wavelet_daubechies, 4);
            work = gsl_wavelet_workspace_alloc(m);

            for (size_t i = 0; i < n; i++) {
                dwtdata[i] = data[i];
            }

            int status = gsl_wavelet_transform_forward(w, dwtdata.data(), 1, m, work);

            if (status != GSL_SUCCESS) {
                printf("Error: wavelets transform failed.\n");
                exit(0);
            }

            for (size_t i = 0; i < n; i++) {
                data[i] = dwtdata[i];
            }

            gsl_wavelet_free(w);
            gsl_wavelet_workspace_free(work);
            

        }


        void preProcess_cdf97(T *data, std::vector<size_t> dims) {
            size_t n=1;
            size_t m_dims=std::array<size_t,3>{1,1,1};
            for (size_t i=0;i<N;i++){
                n*=dims[i];
                m_dims[N-1-i]=dims[i];
            }

            std::vector<double> dwtdata(n, 0);
            for (size_t i = 0; i < n; i++) {
                dwtdata[i] = data[i];
            }

            CDF97 m_cdf;

            m_cdf.take_data(std::move(dwtdata), m_dims);
            auto xforms_xy = num_of_xforms(std::min(m_dims[0], m_dims[1]));
            auto xforms_z = num_of_xforms(m_dims[2]);
            if (xforms_xy == xforms_z)
                m_cdf.dwt3d_dyadic();
            else
                m_cdf.dwt3d_wavelet_packet();


            dwtdata=m_cdf.release_data();

            for (size_t i = 0; i < n; i++) {
                data[i] = dwtdata[i];
            }









        
        

       


        void postProcess(T *data, size_t n) {
            size_t m = n - 1;
            m |= m >> 1;
            m |= m >> 2;
            m |= m >> 4;
            m |= m >> 8;
            m |= m >> 16;
            m++;

            std::vector<double> dwtdata(m, 0);
            gsl_wavelet *w;
            gsl_wavelet_workspace *work;

            w = gsl_wavelet_alloc(gsl_wavelet_daubechies, 4);
            work = gsl_wavelet_workspace_alloc(m);

            for (size_t i = 0; i < n; i++) {
                dwtdata[i] = data[i];
            }

            int status = gsl_wavelet_transform_inverse(w, dwtdata.data(), 1, m, work);

            if (status != GSL_SUCCESS) {
                printf("Error: wavelets transform failed.\n");
                exit(0);
            }

            for (size_t i = 0; i < n; i++) {
                data[i] = dwtdata[i];
            }

            gsl_wavelet_free(w);
            gsl_wavelet_workspace_free(work);

        }

        void postProcess_cdf97(T *data, std::vector<size_t> dims) {
            size_t n=1;
            size_t m_dims=std::array<size_t,3>{1,1,1};
            for (size_t i=0;i<N;i++){
                n*=dims[i];
                m_dims[N-1-i]=dims[i];
            }
            std::vector<double> dwtdata(n, 0);
            for (size_t i = 0; i < n; i++) {
                dwtdata[i] = data[i];
            }

            CDF97 m_cdf;

            m_cdf.take_data(std::move(dwtdata), m_dims);
            auto xforms_xy = num_of_xforms(std::min(m_dims[0], m_dims[1]));
            auto xforms_z = num_of_xforms(m_dims[2]);
            if (xforms_xy == xforms_z)
                m_cdf.idwt3d_dyadic();
            else
                m_cdf.idwt3d_wavelet_packet();


            dwtdata=m_cdf.release_data();

            for (size_t i = 0; i < n; i++) {
                data[i] = dwtdata[i];
            }




        }
        
        
    

       

        
    };
}


#endif
#endif //SZ3_WAVELET_HPP
