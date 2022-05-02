# -*- coding: utf-8 -*-
"""
Created on Sun Dec 12 22:58:54 2021

@author: thinkpad
"""
import time
#demo1
def process_bar(percent, start_str='', end_str='', total_length=0):
    bar = ''.join(["\033[031m%s\033[0m"%'   '] * int(percent * total_length)) + ''
    bar = '\r' + start_str + bar.ljust(total_length) + ' {:0>4.1f}%|'.format(percent*100) + end_str
    print(bar, end='', flush=True)
 
 
for i in range(101):
    time.sleep(0.2)
    end_str = '100%'
    process_bar(i/100, start_str='', end_str=end_str, total_length=15)
