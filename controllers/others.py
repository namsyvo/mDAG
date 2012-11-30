# coding: utf8

def download():
    return response.download(request,db)

def downloads():
    datasets= db().select(db.sampledata.ALL, orderby=db.sampledata.id)
    return dict(datasets=datasets)

def help():
    return dict()

def data():
    return dict()

def references():
    return dict()
