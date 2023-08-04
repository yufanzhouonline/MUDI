
class TAD(object):
    '''
    Class of Topologically Associating Domains (TAD)
    '''
    @staticmethod
    def loadPath(prefix, cellName, res, chrNo):
        '''
        Get the path of TADs.
        '''
        if str(res) == '40000':
            isidsDict = {'ispost':'520001', 'idspost':'240001'}
        elif str(res) == '200000':
            isidsDict = {'ispost':'600001', 'idspost':'400001'}
        elif str(res) == '300000':
            isidsDict = {'ispost':'600001', 'idspost':'600001'}
        elif str(res) == '400000':
            isidsDict = {'ispost':'800001', 'idspost':'800001'}
        elif str(res) == '500000':
            isidsDict = {'ispost':'1000001', 'idspost':'1000001'}
        elif str(res) == '1000000':
            isidsDict = {'ispost':'2000001', 'idspost':'2000001'}
        else:
            isidsDict = {'ispost':'500001', 'idspost':'200001'}
        filePost = '_' + str(res) + '_' + chrno + '.is' + isidsDict['ispost'] + '.ids' + isidsDict['idspost'] + '.insulation.boundaries.bed'
        cellPath = prefix + cellName + '/' + cellName + filePost
        return cellPath
