#!/usr/bin/env node
var assert = require('assert');
var s = require('./build/Release/svd.node');
var A = [
	[1, 2],
	[3, 4],
	[5, 6]
];
console.log('A = %s', JSON.stringify(A));
console.log('---');

// svd using node-svd
var svd = s.svd(A, 0, { U: true, V: false, debug: 1});
console.log('---');
console.log('U = %s', JSON.stringify(svd.U));
console.log('S = %s', JSON.stringify(svd.S));
console.log('V = %s', JSON.stringify(svd.V));

/**
 * Brute-force matrix multiplication
 */
function mult(X, Y){
	assert(X[0].length == Y.length, 'Invalid dimension!');
	// dimensions
	var m = X.length, n = Y[0].length, d = Y.length;
	// multiplication
	var Z = new Array(X.length);
	for(var i = 0; i < X.length; ++i){
		Z[i] = new Array(Y[0].length);
		for(var j = 0; j < Y[0].length; ++j){
			Z[i][j] = X[i].map(function(x, k){
				return x * Y[k][j];
			}).reduce(function(a, b){ return a + b; }, 0);
		}
	}
	return Z;
}

// recomputing A
var S = [[svd.S[0], 0], [0, svd.S[1]]];
var B = mult(svd.U, mult(S, svd.V));
console.log('B = U * diag(S) * Vt\n%s', JSON.stringify(B));
