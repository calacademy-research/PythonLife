#!/usr/bin/env python
 
print 'The top of the greeting_module has been read.'
 
def hello(name):
 greeting = "Hello %s!" % name
 return greeting
 
def ahoy(name):
 greeting = "Ahoy %s!" % name
 return greeting
 
x = 5
 
print 'The bottom of the greeting_module has been read.'
