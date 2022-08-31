extern "C" {
#ifndef SZ3_WAVELET2_HPP
#define SZ3_WAVELET2_HPP

#ifdef ENABLE_GSL 


    #include <gsl/gsl_wavelet.h>




   
      
        void preProcess_t(float *data, size_t n) {
            size_t m = n - 1;
            m |= m >> 1;
            m |= m >> 2;
            m |= m >> 4;
            m |= m >> 8;
            m |= m >> 16;
            m++;

            double * dwtdata=(double *) calloc (m,sizeof(double));
            gsl_wavelet *w;
            gsl_wavelet_workspace *work;

            w = gsl_wavelet_alloc(gsl_wavelet_daubechies, 4);
            work = gsl_wavelet_workspace_alloc(m);

            for (size_t i = 0; i < n; i++) {
                dwtdata[i] = data[i];
            }

            int status = gsl_wavelet_transform_forward(w, dwtdata, 1, m, work);

            if (status != GSL_SUCCESS) {
                printf("Error: wavelets transform failed.\n");
                exit(0);
            }

            for (size_t i = 0; i < n; i++) {
                data[i] = dwtdata[i];
            }
            free(dwtdata);
            gsl_wavelet_free(w);
            gsl_wavelet_workspace_free(work);

        }

        
        
    

 
       


        void postProcess_t(float *data, size_t n) {
            size_t m = n - 1;
            m |= m >> 1;
            m |= m >> 2;
            m |= m >> 4;
            m |= m >> 8;
            m |= m >> 16;
            m++;

            double * dwtdata=(double *) calloc (n,sizeof(double));
            gsl_wavelet *w;
            gsl_wavelet_workspace *work;

            w = gsl_wavelet_alloc(gsl_wavelet_daubechies, 4);
            work = gsl_wavelet_workspace_alloc(m);

            for (size_t i = 0; i < n; i++) {
                dwtdata[i] = data[i];
            }

            int status = gsl_wavelet_transform_inverse(w, dwtdata, 1, m, work);

            if (status != GSL_SUCCESS) {
                printf("Error: wavelets transform failed.\n");
                exit(0);
            }

            for (size_t i = 0; i < n; i++) {
                data[i] = dwtdata[i];
            } 
            free(dwtdata);

            gsl_wavelet_free(w);
            gsl_wavelet_workspace_free(work);

        }
        
        
    

    



#endif
#endif //SZ3_WAVELET2_HPP

    }
