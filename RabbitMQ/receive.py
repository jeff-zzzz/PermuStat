#!/usr/bin/env python
import pika

#connection = pika.BlockingConnection(pika.ConnectionParameters(
#        host='localhost:15672'))
connection = pika.BlockingConnection(pika.ConnectionParameters(
        host='localhost'))
channel = connection.channel()

#channel.queue_declare(queue='com.egi.nolis.86599075-99ee-4a0f-9dda-67e07b84bfd5')
channel.queue_declare(queue='hello')

print ' [*] Waiting for messages. To exit press CTRL+C'

def callback(ch, method, properties, body):
    print " [x] Received %r" % (body,)

#channel.basic_consume(callback,
#                      queue='com.egi.nolis.86599075-99ee-4a0f-9dda-67e07b84bfd5',
#                      no_ack=True)
channel.basic_consume(callback,
                      queue='hello',
                      no_ack=True)

channel.start_consuming()
