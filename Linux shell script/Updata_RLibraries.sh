#!/bin/bash
ssh ng rm -r RFunctions
ssh ng wget -r -np -P RFunctions_tmp "http://viva1109.iptime.org/RFunctions"
ssh ng cp -r ./RFunctions_tmp/viva1109.iptime.org/RFunctions .
ssh ng scp -r ./RFunctions n12:/data/sharedcode/kjkim/
