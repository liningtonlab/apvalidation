import requests
import xml.etree.ElementTree as Et

root_path = r"C:\Users\benle\OneDrive\Documents\Code\ap_validation_production\ap_validation\paramExtract\packages\StructureCreation"


def get_input_data(smiles, input_xml_location = f'{root_path}\input.xml'):
    """
    A function to parse the input xml schema and insert the desired smiles string so that it is ready to
    be sent to the PubChem Standardization tool.
    :param smiles: a smiles string representing a chemical compound
    :return: a xml tree in string format, this tree is what is sent to the PubChem tool
    """
    tree = Et.parse(input_xml_location)
    root = tree.getroot()

    for node in tree.iter():
        if node.tag == 'PCT-Structure_structure_string':
            node.text = smiles

    string_tree = Et.tostring(root, encoding='utf8', method='xml')

    return string_tree


def request_standard(data):
    """
    Make a request to the PubChem standardization tool to standardize the given data.

    :param data: the xml string recieved from get_input_data(smiles)
    :return: a xml tree with a request id sent from PubChem
    """

    response = requests.post('https://pubchem.ncbi.nlm.nih.gov/pug/pug.cgi', data=data)
    response_xml = response.content.decode('utf-8')
    return response_xml


def get_status_data(response_xml, poll_xml_location = f'{root_path}/poll.xml'):
    """
    Use the request id retrieved from request_standard in order to send another request asking
    PubChem for the completed standardization.

    :param response_xml: the xml tree containing the request id
    :return: xml tree to be sent to PubChem. This xml contains the request id
    """
    response_tree = Et.ElementTree(Et.fromstring(response_xml))
    request_id = ''

    for node in response_tree.iter():
        if node.tag == 'PCT-Waiting_reqid':
            request_id = node.text

    input_tree = Et.parse(poll_xml_location)
    root = input_tree.getroot()

    for node in input_tree.iter():
        if node.tag == 'PCT-Request_reqid':
            node.text = request_id

    string_tree = Et.tostring(root, encoding='utf8', method='xml')
    return string_tree


def request_status(status_data):
    """
    Make the second request hoping that PubChem will return the standardized compound.

    :param status_data: The xml tree containing the request id
    :return: xml tree containing the final standard compound.
    """
    response = requests.post('https://pubchem.ncbi.nlm.nih.gov/pug/pug.cgi', data=status_data)
    response_xml = response.content.decode('utf-8')
    return response_xml


def get_standard_struct(smiles):
    """
    A function to call all the above functions in order. One call of this function should
    return a standardized version of the given chemical compound.

    :param smiles: A smiles of the compound to be standardized.
    :return: Standardized smiles.
    """
    input_data = get_input_data(smiles)
    xml_response = request_standard(input_data)
    status_data = get_status_data(xml_response)
    final_xml = request_status(status_data)

    final_tree = Et.ElementTree(Et.fromstring(final_xml))
    final_smiles = "Not found"
    for node in final_tree.iter():
        if node.tag == 'PCT-Structure_structure_string':
            final_smiles = node.text

    return final_smiles



