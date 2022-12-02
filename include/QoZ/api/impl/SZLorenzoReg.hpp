#ifndef SZ3_SZ_LORENZO_REG_HPP
#define SZ3_SZ_LORENZO_REG_HPP

#include "QoZ/compressor/SZGeneralCompressor.hpp"
#include "QoZ/frontend/SZFastFrontend.hpp"
#include "QoZ/frontend/SZGeneralFrontend.hpp"
#include "QoZ/quantizer/IntegerQuantizer.hpp"
#include "QoZ/predictor/ComposedPredictor.hpp"
#include "QoZ/predictor/LorenzoPredictor.hpp"
#include "QoZ/predictor/RegressionPredictor.hpp"
#include "QoZ/predictor/PolyRegressionPredictor.hpp"
#include "QoZ/predictor/ZeroPredictor.hpp"
#include "QoZ/lossless/Lossless_zstd.hpp"
#include "QoZ/utils/Iterator.hpp"
#include "QoZ/utils/Statistic.hpp"
#include "QoZ/utils/Extraction.hpp"
#include "QoZ/utils/QuantOptimizatioin.hpp"
#include "QoZ/utils/Config.hpp"
#include "QoZ/def.hpp"
#include <cmath>
#include <cstdlib>
#include <memory>


template<class T, QoZ::uint N, class Quantizer, class Encoder, class Lossless>
std::shared_ptr<QoZ::concepts::CompressorInterface<T>>
make_lorenzo_regression_compressor(const QoZ::Config &conf, Quantizer quantizer, Encoder encoder, Lossless lossless) {
    std::vector<std::shared_ptr<QoZ::concepts::PredictorInterface<T, N>>> predictors;

    int methodCnt = (conf.lorenzo + conf.lorenzo2 + conf.regression + conf.regression2);
    int use_single_predictor = (methodCnt == 1);
    if (methodCnt == 0) {
        printf("All lorenzo and regression methods are disabled.\n");
        exit(0);
    }
    if (conf.lorenzo) {
        std::vector<double> coeffs;
        if(conf.useCoeff)
            coeffs=conf.lorenzo1_coeffs;
        if (use_single_predictor) {
            return QoZ::make_sz_general_compressor<T, N>(
                    QoZ::make_sz_general_frontend<T, N>(conf, QoZ::LorenzoPredictor<T, N, 1>(conf.absErrorBound,coeffs), quantizer),
                    encoder, lossless);
        } else {
            predictors.push_back(std::make_shared<QoZ::LorenzoPredictor<T, N, 1>>(conf.absErrorBound,coeffs));
        }
    }
    if (conf.lorenzo2) {
        std::vector<double> coeffs;
        if(conf.useCoeff)
            coeffs=conf.lorenzo2_coeffs;
        if (use_single_predictor) {
            return QoZ::make_sz_general_compressor<T, N>(
                    QoZ::make_sz_general_frontend<T, N>(conf, QoZ::LorenzoPredictor<T, N, 2>(conf.absErrorBound,coeffs), quantizer),
                    encoder, lossless);
        } else {
            predictors.push_back(std::make_shared<QoZ::LorenzoPredictor<T, N, 2>>(conf.absErrorBound,coeffs));
        }
    }
    if (conf.regression) {
        if (use_single_predictor) {
            return QoZ::make_sz_general_compressor<T, N>(
                    QoZ::make_sz_general_frontend<T, N>(conf, QoZ::RegressionPredictor<T, N>(conf.blockSize, conf.absErrorBound),
                                                       quantizer), encoder, lossless);
        } else {
            predictors.push_back(std::make_shared<QoZ::RegressionPredictor<T, N>>(conf.blockSize, conf.absErrorBound));
        }
    }

    if (conf.regression2) {
        if (use_single_predictor) {
            return QoZ::make_sz_general_compressor<T, N>(
                    QoZ::make_sz_general_frontend<T, N>(conf, QoZ::PolyRegressionPredictor<T, N>(conf.blockSize, conf.absErrorBound),
                                                       quantizer), encoder, lossless);
        } else {
            predictors.push_back(std::make_shared<QoZ::PolyRegressionPredictor<T, N>>(conf.blockSize, conf.absErrorBound));
        }
    }
    return QoZ::make_sz_general_compressor<T, N>(
            QoZ::make_sz_general_frontend<T, N>(conf, QoZ::ComposedPredictor<T, N>(predictors),
                                               quantizer), encoder, lossless);
}


template<class T, QoZ::uint N>
char *SZ_compress_LorenzoReg(QoZ::Config &conf, T *data, size_t &outSize) {

    assert(N == conf.N);
    assert(conf.cmprAlgo == QoZ::ALGO_LORENZO_REG);
    //QoZ::calAbsErrorBound(conf, data);

    char *cmpData;
    auto quantizer = QoZ::LinearQuantizer<T>(conf.absErrorBound, conf.quantbinCnt / 2);
    if (N == 3 and !conf.regression2 and !conf.useCoeff) {
        // use fast version for 3D
        auto sz = QoZ::make_sz_general_compressor<T, N>(QoZ::make_sz_fast_frontend<T, N>(conf, quantizer), QoZ::HuffmanEncoder<int>(),
                                                       QoZ::Lossless_zstd());
        cmpData = (char *) sz->compress(conf, data, outSize);
    } else {
        auto sz = make_lorenzo_regression_compressor<T, N>(conf, quantizer, QoZ::HuffmanEncoder<int>(), QoZ::Lossless_zstd());
        //std::cout<<"lor1"<<std::endl;
        cmpData = (char *) sz->compress(conf, data, outSize);
    }
    return cmpData;
}


template<class T, QoZ::uint N>
void SZ_decompress_LorenzoReg(const QoZ::Config &theconf, char *cmpData, size_t cmpSize, T *decData) {
    QoZ::Config conf(theconf);
    assert(conf.cmprAlgo == QoZ::ALGO_LORENZO_REG);
    QoZ::uchar const *cmpDataPos = (QoZ::uchar *) cmpData;
    QoZ::LinearQuantizer<T> quantizer;
    if(!conf.wavelet){

        
        if (N == 3 and !conf.regression2 and !conf.useCoeff) {
            // use fast version for 3D
            auto sz = QoZ::make_sz_general_compressor<T, N>(QoZ::make_sz_fast_frontend<T, N>(conf, quantizer),
                                                           QoZ::HuffmanEncoder<int>(), QoZ::Lossless_zstd());
            sz->decompress(cmpDataPos, cmpSize, decData);
           
        } else {
            auto sz = make_lorenzo_regression_compressor<T, N>(conf, quantizer, QoZ::HuffmanEncoder<int>(), QoZ::Lossless_zstd());
            sz->decompress(cmpDataPos, cmpSize, decData);
           
        }
    }
    else{
        std::vector<size_t> ori_dims=conf.dims;
        size_t ori_num=conf.num;
        if(conf.external_wave){
            conf.setDims(conf.coeffs_dims.begin(),conf.coeffs_dims.end());
            
        }

        size_t first =conf.firstSize;
        size_t second=cmpSize-conf.firstSize;
        if (N == 3 and !conf.regression2 and !conf.useCoeff) {
            // use fast version for 3D
            auto sz = QoZ::make_sz_general_compressor<T, N>(QoZ::make_sz_fast_frontend<T, N>(conf, quantizer),
                                                           QoZ::HuffmanEncoder<int>(), QoZ::Lossless_zstd());
            sz->decompress(cmpDataPos, first, decData);
           
        } else {
            auto sz = make_lorenzo_regression_compressor<T, N>(conf, quantizer, QoZ::HuffmanEncoder<int>(), QoZ::Lossless_zstd());
            sz->decompress(cmpDataPos, first, decData);
           
        }

        /*
        if(conf.transformation==1){
            for(size_t i=0;i<conf.num;i++)
                decData[i]=QoZ::logit<double>(decData[i]);
        }
        else if(conf.transformation==2){
            for(size_t i=0;i<conf.num;i++)
                decData[i]=QoZ::arctanh<double>(decData[i]);
        } 
        */



         //QoZ::writefile<T>("waved.qoz.dec.logit", decData, conf.num);
        if(conf.external_wave){

            char s1[100]=itoa(conf.pid);
            char s2[]="_external_dec_wave_coeffs_dec.tmp";
            strcat(s1,s2);
            QoZ::writefile(s1, decData, conf.num);

            char command[120] = "python coeff_idwt.py ";//still need slice.pkl wave_type.txt wave_size.dat, or pickle all metadata into one file.
            strcat(command,s1);
            system(command);

           conf.setDims(ori_dims.begin(),ori_dims.end());
          
            delete []decData;
            decData=new T[conf.num];
            char s3[100]=itoa(conf.pid);
            char s4[]="_external_deccoeff_idwt.tmp";
            strcat(s3,s4);
            QoZ::readfile<T>(s3, conf.num, decData);
        }

        else{
            QoZ::Wavelet<T,N> wlt;

            wlt.postProcess_cdf97(decData,conf.dims);
        }
        //QoZ::writefile<T>("waved.qoz.dec.idwt", decData, conf.num);


        
       
       
        T *offsets =new T [conf.num];
        

        QoZ::Config newconf(conf.num);

        //newconf.blockSize=32768;

        if (conf.offsetPredictor ==0){
            auto quantizer = QoZ::LinearQuantizer<T>(newconf.absErrorBound, newconf.quantbinCnt / 2);
            auto sz2 = QoZ::make_sz_general_compressor<T, 1>(QoZ::make_sz_general_frontend<T, 1>(newconf, QoZ::ZeroPredictor<T, 1>(), quantizer), QoZ::HuffmanEncoder<int>(),
                                                                       QoZ::Lossless_zstd());
           
            
       
             sz2->decompress(cmpDataPos+first,second,offsets);
        }

        else if (conf.offsetPredictor ==1){
            newconf.lorenzo = true;
            newconf.lorenzo2 = true;
            newconf.regression = false;
            newconf.regression2 = false;
            newconf.openmp = false;
            newconf.blockSize = 16;//original 5
            newconf.quantbinCnt = 65536 * 2;

            auto quantizer = QoZ::LinearQuantizer<T>(newconf.absErrorBound, newconf.quantbinCnt / 2);
            auto sz2 = make_lorenzo_regression_compressor<T, 1>(newconf, quantizer, QoZ::HuffmanEncoder<int>(), QoZ::Lossless_zstd());
           
            
       
              sz2->decompress(cmpDataPos+first,second,offsets);
        }
        else if (conf.offsetPredictor == 2){
            newconf.setDims(conf.dims.begin(),conf.dims.end());
            newconf.lorenzo = true;
            newconf.lorenzo2 = true;
            newconf.regression = false;
            newconf.regression2 = false;
            newconf.openmp = false;
            newconf.blockSize = 5;
            newconf.quantbinCnt = 65536 * 2;

            auto quantizer = QoZ::LinearQuantizer<T>(newconf.absErrorBound, newconf.quantbinCnt / 2);
            auto sz2 = make_lorenzo_regression_compressor<T, N>(newconf, quantizer, QoZ::HuffmanEncoder<int>(), QoZ::Lossless_zstd());
           
            
       
            sz2->decompress(cmpDataPos+first,second,offsets);
        }

        else if (conf.offsetPredictor == 3){
            newconf.interpAlgo=QoZ::INTERP_ALGO_CUBIC;
            newconf.interpDirection=0;

           
            auto sz2 = QoZ::SZInterpolationCompressor<T, 1, QoZ::LinearQuantizer<T>, QoZ::HuffmanEncoder<int>, QoZ::Lossless_zstd>(
            QoZ::LinearQuantizer<T>(newconf.absErrorBound),
            QoZ::HuffmanEncoder<int>(),
            QoZ::Lossless_zstd());
            
       
            sz2.decompress(cmpDataPos+first,second,offsets);
        }

        else if (conf.offsetPredictor == 4){
            
            newconf.setDims(conf.dims.begin(),conf.dims.end());
            newconf.interpAlgo=QoZ::INTERP_ALGO_CUBIC;
            newconf.interpDirection=0;
           
            auto sz2 = QoZ::SZInterpolationCompressor<T, N, QoZ::LinearQuantizer<T>, QoZ::HuffmanEncoder<int>, QoZ::Lossless_zstd>(
            QoZ::LinearQuantizer<T>(newconf.absErrorBound),
            QoZ::HuffmanEncoder<int>(),
            QoZ::Lossless_zstd());
            
       
            sz2.decompress(cmpDataPos+first,second,offsets);
        }
        //QoZ::writefile<T>("waved.qoz.dec.offset", offsets, conf.num);


        
       

        for(size_t i=0;i<conf.num;i++)
            decData[i]+=offsets[i];
      

        //delete [] cmpDataFirst;
        //delete [] cmpDataSecond;
        delete [] offsets;

        
    }




}

#endif