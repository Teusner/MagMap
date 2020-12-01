# MagMap

[![forthebadge](https://forthebadge.com/images/badges/made-with-python.svg)](https://forthebadge.com)
[![forthebadge](https://forthebadge.com/images/badges/60-percent-of-the-time-works-every-time.svg)](https://forthebadge.com)
[![forthebadge](https://forthebadge.com/images/badges/powered-by-coffee.svg)](https://forthebadge.com)

## Introduction

This project is a student project about the mapping of a field using a magnetometer. This magnetometer need to be trailed on a sledge behind the robot in order to have better measurements because the trailing robot generate magnetic perturbations.

## Some already made content

Using feedback linearization, we are able to generate the control vector ![](https://render.githubusercontent.com/render/math?math=\begin{bmatrix}u_1\\u_2\end{bmatrix}) to control the speed and the bearing of the robot. This robot is then able to follow any wanted trajectory as the following lissajous curve :

![](./documentation/imgs/simulator_lissajous.png)