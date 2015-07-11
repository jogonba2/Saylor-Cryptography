#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  my_crypto.py
#  
#  Copyright 2015 Overxflow13

from math import sqrt,floor,log
from random import randint

def exp(a,b):
	""" Pow using binary exponentiation """
	if b==1: return a
	else:
		if a%2==0: return exp(a,b/2)*exp(a,b/2)
		else: return a * exp(a,b-1)

def primality_test_6k(n):
	""" Test primality using 6k+-i method """
	c = 1
	if n%2==0 or n%3==0: return False
	else:
		for k in xrange(1,int(sqrt(n))):
			p = 6*k
			if n%(p+1)==0 or n%(p-1)==0: c += 1
	if c==1: return True
	else: return False

def mcd(a,b):
	""" Computes mcd using extended iterated euclides algorithm """
	s,t,p,m = 1,0,0,1
	while b!=0:
		q,r = a/b,a%b
		a,s,t,b,p,m = b,p,m,r,s-p*q,t-m*q
	return a

def coprime(n):
	""" Finds coprime to n +- naive iterated algorithm """
	for i in xrange(2,int(sqrt(n))):
		if mcd(n,i)==1: return i
	return n
	
def fermat_primality_test(n):
	""" Computes fermat primality test """
	if (coprime(n)**(n-1))%n != 1: return False
	return True

def pollards_rho(n,it,f):
	""" Computes pollards-rh0 factorization """
	x,y,d,i = randint(1,n-3),randint(1,n-1),1,0
	while d==1 and i<it:
		x,y = f(x),f(f(y))
		d   = mcd(abs(x-y),n)
		i   += it
	return d

#Thanks to http://martin-thoma.com/calculate-legendre-symbol/ #
def calculateLegendre(a, p):
    """ Calculate the legendre symbol (a, p) with p is prime. """
    if a >= p or a < 0:
        return calculateLegendre(a % p, p)
    elif a == 0 or a == 1:
        return a
    elif a == 2:
        if p%8 == 1 or p%8 == 7:
            return 1
        else:
            return -1
    elif a == p-1:
        if p%4 == 1:
            return 1
        else:
            return -1
    elif not primality_test_6k(a):
        factors = factorize(a)
        product = 1
        for pi in factors:
            product *= calculateLegendre(pi, p)
        return product
    else:
        if ((p-1)/2)%2==0 or ((a-1)/2)%2==0:
            return calculateLegendre(p, a)
        else:
            return (-1)*calculateLegendre(p, a)

# Thanks to Alex Martelli #
def perfect_square(apositiveint):
  x = apositiveint // 2
  seen = set([x])
  while x * x != apositiveint:
    x = (x + (apositiveint // x)) // 2
    if x in seen: return False
    seen.add(x)
  return True
  
def shanks(n,k):
	""" Computes Shanks factorization """
	# First Phase #
	p0,q0 = floor(sqrt(k*n)),1
	q1 = (k*n)-(p0**2)
	actQ,antQ,p = q1,q0,p0
	while perfect_square(actQ)==False:
		b    = floor((floor(sqrt(k*n))+p)/actQ)
		antP,p = p,(b*actQ)-p
		antQ,actQ = actQ,antQ + b*(antP-p)
		#print p,antQ,actQ
	# Second Phase #
	b      = floor((floor(sqrt(k*n))-p)/(sqrt(actQ)))
	p0     = (b*sqrt(actQ))+p
	antQ   = sqrt(actQ)
	actQ   = ((k*n)-(p0**2))/antQ
	p,antP = p0,-1
	print p0,antQ,actQ
	while antP!=p:
		b    = floor((floor(sqrt(k*n))+p)/actQ)
		antP,p = p,(b*actQ)-p
		antQ,actQ = actQ,antQ + b*(antP-p)
		#print p,antQ,actQ
	m = mcd(n,p)
	if m!=1 and m!=n: return False
	return True


def breakdown_number(n):
	""" Naive number breakdown """
	for d in xrange(n):
		for s in xrange(1,int(log(n,2)+1)):
			if d*pow(2,s)+1==n: return (d,s)
			
	
def strong_pseudoprime(n,a):
	v = breakdown_number(n)
	if v == None: return False
	else:
		d,s = v[0],v[1]	
		if (a**d)%n == 1%n: return True
		else:
			for r in xrange(s-1):
				if (a**(d*(2**r)))%n == -1%n: return True
	return False
	
def solovay_strassen(n):
	""" Computes solovay-strassen test (using legendre not jacobi at this moment) """
	a = randint(1,n-1)
	if mcd(a,n)!=1: return False
	else:
		if calculateLegendre(a,n)%n != (a**((n-1)/2))%n: return False
	return True
	
def miller_rabin(n):
	""" Computes miller-rabin algorithm """
	a = randint(2,n-1)
	if mcd(a,n)!=1: return False
	if mcd(a,n)==1:
		# False -> a^((2^s)*n)%n == -1%n with s € {0,...,t-1} / n = (2^t)*p #
		if (a**n)%n!=1%n or (a**n)%n!=-1%n or False: return False
	return True
	
def criba_eratostenes(n):
	l,mark = [],[]
	for i in xrange(2,int(sqrt(n))+1):
		if i not in mark:
			for j in xrange(i,int(n/i)+1): mark.append(i*j)
	return [x for x in xrange(2,n+1) if x not in mark]
