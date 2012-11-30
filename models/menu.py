response.title = 'mDAG'
response.subtitle = T('analysis of multiple-treatment microarray data')

response.menu = [(T('Help'), False, URL(request.application,'others','help'),
                  [
                      (T('Data Format'), False, URL(request.application,'others','data')),
                      (T('Downloads'), False, URL(request.application,'others','downloads')),
                      (T('References'), False, URL(request.application,'others','references')),
                  ]),
                 (T('Data'), False, URL(request.application,'default','data')),
                 (T('Analyze'), False, URL(request.application,'default','analyze')),
                 (T('Results'), False, URL(request.application,'default','results'))
                 ]
if 'auth' in globals() and auth.user:
    if auth.user.id == 1:
        response.menu.append((T('Administration'), False, URL(request.application,'appadmin','index')))
