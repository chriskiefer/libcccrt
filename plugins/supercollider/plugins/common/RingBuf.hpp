/*
Part of the Fluid Corpus Manipulation Project (http://www.flucoma.org/)
Copyright 2017-2019 University of Huddersfield.
Licensed under the BSD-3 License.
See license.md file in the project root for full license information.
This project has received funding from the European Research Council (ERC)
under the European Union’s Horizon 2020 research and innovation programme
(grant agreement No 725899).

*/

#pragma once

#include <Eigen/Core>
#include <cmath>



template <typename T>
class RingBuf {
public:
    using ringBufType = Eigen::Array<T,Eigen::Dynamic,1>;
    using winBufType = Eigen::Array<T,Eigen::Dynamic,1>;

    RingBuf() {
        buf.resize(512);
        buf.fill(0);
        currBuf.resize(512);
    }

    void setSize(size_t size) {
        buf.resize(size);
        buf.fill(0);
        idx=0;
    }

    void push(T x) {
        buf[idx] = x;
        idx++;
        if (idx==buf.size()) {
            idx=0;
        }
    }
    
    size_t size() {return buf.size();}
    
    winBufType& getBuffer(const size_t winSize, const size_t offset = 0) {
        if (winSize != buf.size())
          currBuf.resize(winSize);
        
        int targidx=0;
        const size_t winAndOffset = winSize + offset;
        if (idx > winAndOffset) {
            for(size_t i=idx-winAndOffset; i < idx-offset; i++, targidx++) {
                currBuf[targidx] = buf[i];
            }
        }else{
            size_t winStart = buf.size()-(winAndOffset-idx);
            size_t winEnd1 = winStart + winSize;
            size_t winEnd2 = 0;
            if (winEnd1 > buf.size()) {
                winEnd2 = winEnd1 - buf.size();
                winEnd1 = buf.size();
            }
            //first chunk
            for(size_t i=winStart; i < winEnd1; i++, targidx++) {
                currBuf[targidx] = buf[i];
            }
            //second chunk
            for(size_t i=0; i < winEnd2; i++, targidx++) {
                currBuf[targidx] = buf[i];
            }
        }
        return currBuf;
    }
    
private:
    ringBufType buf;
    size_t idx=0;
    winBufType currBuf;

};

