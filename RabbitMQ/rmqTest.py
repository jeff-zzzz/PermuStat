#!/usr/bin/python

from RabbitMQ import Connection
import numpy

c = Connection.Connection('localhost')
c.connect()

ar = numpy.arange(2100.0) / 2100.0
ar2 = numpy.arange(2100.0*3.0) / (2100.0*3.0)
ar3 = numpy.arange(256.0) / 128.0 - 1.0

c.sendOrientedData(ar, "Oriented Test", "123")
c.sendTriplesData(ar2, "Triples Test", "1234")

c.sendEEGData(ar3, "Oriented Test", "123")
c.sendEEGData(ar3, "Triples Test", "1234")
c.disconnect()
