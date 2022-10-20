import cccrt

#ETC 

a = array([random.randint(0,2) for x in arange(100)])
print(a)
print(cccrt.ETC(a))

#Joint ETC

a = arange(100) % 20
print(a)
b = array(sin(a) * 10, dtype=int)
print(b)
print(cccrt.ETCJoint(a,b))

#dynamic CC

a = arange(100) % 20
a = a + array([random.randint(0,4) for x in arange(100)])
print(a)
print(cccrt.dynamicCC(a, 10, 10, 1))

#dynamic joint CC

a = arange(100) % 20
print(a)
b = array(sin(a) * 10, dtype=int)
print(b)
print(cccrt.dynamicCCJoint(a,b, 5,5,5))

###  CCC

a = arange(100) % 3
a = a + array([random.randint(0,10) for x in arange(100)])
print(a)
b = array(sin(a) * 10, dtype=int)
print(b)
print(cccrt.CCCausality(a,b, 5,5,5))


# Lempel Ziv

a = arange(10) % 4
print(a)
print(cccrt.LZ(a))  #un-normalised
print(cccrt.NLZ(a)) #normalised


