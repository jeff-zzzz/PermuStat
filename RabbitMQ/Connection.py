import pika
#import numpy

class Connection:
    def __init__(self, server):
        self.connection = None
        self.channel = None
        self.exchange = 'com.egi.nolis.geosource.exchange'
        self.queue = 'PythonRMQ'
        self.params = pika.ConnectionParameters(server)
        
    def connect(self):
        self.connection = pika.BlockingConnection(self.params)
        self.channel = self.connection.channel()
        self.channel.exchange_declare(exchange=self.exchange, type='fanout', durable=True)

    def disconnect(self):
        self.connection.close()
        self.channel = None
        
    def sendGenericMessage(self, msg):
        self.channel.basic_publish(exchange=self.exchange, routing_key='', body=msg, properties=pika.BasicProperties(content_type='text/plain', delivery_mode=1))

    def sendTriplesData(self, dataAr, conditionName, timeString):
        msg = "<NolisCommandSet>\n"
        msg += "<NolisSetTriplesData NumValues=\"" + str(dataAr.size) + "\" Encoding=\"ascii\" Condition=\"" + conditionName + "\" Time=\"" + timeString + "\" DataType=\"float\">\n"
        msg += "  <Data>\n"
        for i in range(dataAr.flatten().shape[0]):
            msg += str(dataAr.flatten()[i]) + " "
        msg += "\n  </Data>\n"
        msg += "</NolisSetTriplesData>\n"
        msg += "</NolisCommandSet>"
        self.sendGenericMessage(msg)

    def sendOrientedData(self, dataAr, conditionName, timeString):
        msg = "<NolisCommandSet>\n"
        msg += "  <NolisSetDipoleData NumValues=\"" + str(dataAr.size) + "\" Encoding=\"ascii\" Condition=\"" + conditionName + "\" Time=\"" + timeString + "\" DataType=\"float\" MeshID=\"-1\">\n"
        msg += "    <Data>\n"
        for i in range(dataAr.flatten().shape[0]):
            msg += str(dataAr.flatten()[i]) + " "
        msg += "\n    </Data>\n"
        msg += "  </NolisSetDipoleData>\n"
        msg += "</NolisCommandSet>"
        self.sendGenericMessage(msg)

    # Send over just data that belongs on EEG sensors
    def sendEEGData(self, dataAr, conditionName, timeString):
        msg = "<NolisCommandSet>\n"
        msg += "<NolisEEGData NumValues=\"" + str(dataAr.size) + "\" Encoding=\"ascii\" Condition=\"" + conditionName + "\" Time=\"" + timeString + "\" DataType=\"float\">\n"
        msg += "  <Data>\n"
        for i in range(dataAr.flatten().shape[0]):
            msg += str(dataAr.flatten()[i]) + " "
        msg += "\n  </Data>\n"
        msg += "</NolisEEGData>\n"
        msg += "</NolisCommandSet>"
        self.sendGenericMessage(msg)
