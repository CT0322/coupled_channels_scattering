/*
 * Momentum_lookup.h
 *
 *  Created on: Mar 26, 2013
 *      Author: knippsch
 */

#ifndef MOMENTUMLOOKUP_H_
#define MOMENTUMLOOKUP_H_

#include <iostream>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <typeinfo>

#include <Eigen/Dense>
#include <Eigen/Core>
#include <Eigen/SparseCore>



void Momentum_lookup(Eigen::Vector3i, int, int, Eigen::Vector2i**&, int&);

#endif /* MOMENTUMLOOKUP_H_ */
