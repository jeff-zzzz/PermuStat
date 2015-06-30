
def ns2pyTime(nsTime):
    from datetime import datetime
    nsDate = nsTime[0:10]
    nsTime0 = nsTime[11:26]
    nsTime00 = nsDate + " " + nsTime0
    pyTime = datetime.strptime(nsTime00, '%Y-%m-%d %H:%M:%S.%f')
    return pyTime

## beginTime0 = ns2pyTime(beginTime)
## eventTime0 = ns2pyTime(eventTime)
## (eventTime0-beginTime0).microseconds

