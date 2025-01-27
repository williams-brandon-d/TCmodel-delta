# -*- coding: utf-8 -*-
"""
Created on Mon Nov 29 13:28:29 2021

@author: brndn
"""

def vecToNumpy(trec,vrec_tc0,vrec_tc1,vrec_re0,vrec_re1):
    time = trec.as_numpy()
    v_tc0 = vrec_tc0.as_numpy()
    v_tc1 = vrec_tc1.as_numpy()
    v_re0 = vrec_re0.as_numpy()
    v_re1 = vrec_re1.as_numpy()
    return time,v_tc0,v_tc1,v_re0,v_re1
