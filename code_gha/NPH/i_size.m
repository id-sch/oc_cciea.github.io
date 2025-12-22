%%% i_size.m --- 
%% 
%% Filename: i_size.m
%% Description: 
%% Author: Isaac  Schroeder
%% Maintainer: 
%% Created: Thu Jan 10 13:53:41 2008
%% Version: $Id$
%% Last-Updated: Thu Jan 10 13:54:51 2008
%%           By: Isaac  Schroeder
%%     Update #: 1
%% Keywords: 
%% Compatibility: 
%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 
%%% Commentary: size of a 1d matrix
%% 
%% 
%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 
%%% Change log:
%% 
%% RCS $Log$
%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 
%%% Code:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function nop=i_size(data)
    
    
nop = prod(size(data));