import xml.etree.ElementTree as ET
import requests
import gzip
from collections import defaultdict
import json

class MeSHHierarchyExtractor:
    def __init__(self):
        self.mesh_hierarchy = defaultdict(list)
        self.descriptor_to_tree = defaultdict(list)
        self.tree_to_descriptor = {}
        self.parent_child_relations = defaultdict(set)
        self.descriptor_details = {}  # Store complete MeSH information
        self.ui_to_tree = defaultdict(list)  # MeSH ID to tree numbers
        self.tree_to_ui = {}  # Tree number to MeSH ID
        
    def download_mesh_xml(self, year=2020):
        """
        Download MeSH XML files for the specified year
        """
        base_url = f"https://nlmpubs.nlm.nih.gov/projects/mesh/MESH_FILES/xmlmesh/"
        
        # For 2020, the descriptor file would be desc2020.xml
        descriptor_url = f"{base_url}desc{year}.xml"
        
        print(f"Downloading MeSH descriptors for {year}...")
        try:
            response = requests.get(descriptor_url)
            response.raise_for_status()
            return response.content
        except requests.exceptions.RequestException as e:
            print(f"Error downloading: {e}")
            print("Try downloading manually from: https://nlmpubs.nlm.nih.gov/projects/mesh/")
            return None
    
    def parse_mesh_xml(self, xml_content):
        """
        Parse MeSH XML and extract hierarchical relationships with key identifiers
        """
        if isinstance(xml_content, bytes):
            xml_content = xml_content.decode('utf-8')
            
        root = ET.fromstring(xml_content)
        
        # Store complete descriptor information
        self.descriptor_details = {}
        
        for descriptor_record in root.findall('.//DescriptorRecord'):
            # Get MeSH Unique ID (most important - primary identifier)
            descriptor_ui_elem = descriptor_record.find('DescriptorUI')
            if descriptor_ui_elem is None:
                continue
            descriptor_ui = descriptor_ui_elem.text
            
            # Get descriptor name
            descriptor_name_elem = descriptor_record.find('DescriptorName/String')
            if descriptor_name_elem is None:
                continue
            descriptor_name = descriptor_name_elem.text
            
            # Get annotation/definition (scope note)
            annotation = ""
            annotation_elem = descriptor_record.find('Annotation')
            if annotation_elem is not None:
                annotation = annotation_elem.text
            
            # Alternative: Get concept scope note
            if not annotation:
                concept_elem = descriptor_record.find('.//Concept/ScopeNote')
                if concept_elem is not None:
                    annotation = concept_elem.text
            
            # Store complete descriptor details
            self.descriptor_details[descriptor_ui] = {
                'mesh_id': descriptor_ui,
                'name': descriptor_name,
                'definition': annotation or "",
                'tree_numbers': []
            }
            
            # Get tree numbers (hierarchy positions)
            tree_number_list = descriptor_record.find('TreeNumberList')
            if tree_number_list is not None:
                for tree_number_elem in tree_number_list.findall('TreeNumber'):
                    tree_number = tree_number_elem.text
                    
                    # Add to descriptor details
                    self.descriptor_details[descriptor_ui]['tree_numbers'].append(tree_number)
                    
                    # Store mappings (keeping original functionality)
                    self.descriptor_to_tree[descriptor_name].append(tree_number)
                    self.tree_to_descriptor[tree_number] = descriptor_name
                    
                    # Store UI to tree mapping
                    if not hasattr(self, 'ui_to_tree'):
                        self.ui_to_tree = defaultdict(list)
                    self.ui_to_tree[descriptor_ui].append(tree_number)
                    
                    # Store tree to UI mapping
                    if not hasattr(self, 'tree_to_ui'):
                        self.tree_to_ui = {}
                    self.tree_to_ui[tree_number] = descriptor_ui
                    
                    # Extract parent-child relationships
                    self._extract_hierarchy_relations(tree_number, descriptor_name)
    
    def _extract_hierarchy_relations(self, tree_number, descriptor_name):
        """
        Extract parent-child relationships from tree numbers
        Tree numbers like C01.123.456 indicate hierarchy levels
        """
        parts = tree_number.split('.')
        
        # Build parent tree numbers
        for i in range(1, len(parts)):
            parent_tree = '.'.join(parts[:i])
            child_tree = '.'.join(parts[:i+1])
            
            if parent_tree in self.tree_to_descriptor and child_tree in self.tree_to_descriptor:
                parent_name = self.tree_to_descriptor[parent_tree]
                child_name = self.tree_to_descriptor[child_tree]
                self.parent_child_relations[parent_name].add(child_name)
    
    def get_label_hierarchy_for_micol(self):
        """
        Generate label-label similarity pairs for MICoL enhancement
        Returns pairs of (label1, label2) that are hierarchically related
        """
        similar_pairs = []
        
        # Parent-child relationships (strong similarity)
        for parent, children in self.parent_child_relations.items():
            for child in children:
                similar_pairs.append((parent, child, "parent_child"))
        
        # Sibling relationships (medium similarity)
        for parent, children in self.parent_child_relations.items():
            children_list = list(children)
            for i, child1 in enumerate(children_list):
                for child2 in children_list[i+1:]:
                    similar_pairs.append((child1, child2, "siblings"))
        
    def get_mesh_by_id(self, mesh_id):
        """
        Get complete MeSH information by Unique ID (like D016667)
        Returns the highlighted key-value pairs from the browser
        """
        if mesh_id in self.descriptor_details:
            details = self.descriptor_details[mesh_id]
            return {
                'MeSH_Unique_ID': mesh_id,
                'Tree_Numbers': details['tree_numbers'],
                'Term_Name': details['name'],
                'Definition': details['definition']
            }
        return None
    
    def search_mesh_by_name(self, term_name):
        """
        Find MeSH ID by term name (like "Bacterial Capsules")
        """
        for mesh_id, details in self.descriptor_details.items():
            if details['name'].lower() == term_name.lower():
                return self.get_mesh_by_id(mesh_id)
        return None
    
    def get_hierarchy_relations_by_id(self, mesh_id):
        """
        Get parent-child relationships using MeSH Unique ID
        """
        if mesh_id not in self.descriptor_details:
            return None
            
        term_name = self.descriptor_details[mesh_id]['name']
        tree_numbers = self.descriptor_details[mesh_id]['tree_numbers']
        
        parents = []
        children = list(self.parent_child_relations.get(term_name, []))
        
        # Find parents by looking at shorter tree numbers
        for tree_num in tree_numbers:
            parts = tree_num.split('.')
            if len(parts) > 1:
                parent_tree = '.'.join(parts[:-1])
                if parent_tree in self.tree_to_ui:
                    parent_id = self.tree_to_ui[parent_tree]
                    parent_name = self.descriptor_details[parent_id]['name']
                    parents.append({'id': parent_id, 'name': parent_name})
        
        return {
            'mesh_id': mesh_id,
            'term_name': term_name,
            'parents': parents,
            'children': [{'name': child} for child in children],
            'tree_numbers': tree_numbers
        }
    
    def save_hierarchy_data(self, output_file="mesh_hierarchy_v2.json"):
        """
        Save extracted hierarchy data to JSON file with complete MeSH information
        """
        hierarchy_data = {
            "descriptor_details": self.descriptor_details,  # Complete MeSH info with IDs
            "descriptor_to_tree": dict(self.descriptor_to_tree),
            "tree_to_descriptor": self.tree_to_descriptor,
            "ui_to_tree": dict(self.ui_to_tree),  # MeSH ID to tree mapping
            "tree_to_ui": self.tree_to_ui,  # Tree to MeSH ID mapping
            "parent_child_relations": {k: list(v) for k, v in self.parent_child_relations.items()},
            "similar_pairs": self.get_label_hierarchy_for_micol()
        }
        
        with open(output_file, 'w') as f:
            json.dump(hierarchy_data, f, indent=2)
        
        print(f"Hierarchy data saved to {output_file}")
        print(f"Sample MeSH entries with IDs:")
        
        # Show sample entries with the key fields you highlighted
        for i, (mesh_id, details) in enumerate(self.descriptor_details.items()):
            if i >= 3:  # Show first 3 examples
                break
            print(f"  MeSH ID: {mesh_id}")
            print(f"  Name: {details['name']}")
            print(f"  Tree Numbers: {details['tree_numbers']}")
            if details['definition']:
                print(f"  Definition: {details['definition'][:100]}...")
            print()
        
        return hierarchy_data

def main():
    # Example usage
    extractor = MeSHHierarchyExtractor()
    
    try:
        with open('desc2025.xml', 'r', encoding='utf-8') as f:
            xml_content = f.read()
        
        print("Parsing MeSH XML...")
        extractor.parse_mesh_xml(xml_content)
        
        print("Extracting hierarchy relationships...")
        hierarchy_data = extractor.save_hierarchy_data()
        
        # Print some statistics
        print(f"\nStatistics:")
        print(f"Total descriptors: {len(extractor.tree_to_descriptor)}")
        print(f"Total parent-child relations: {sum(len(children) for children in extractor.parent_child_relations.values())}")
        # ßprint(f"Total similar pairs for training: {len(hierarchy_data['similar_pairs'])}")
        
        # Show example relationships
        print(f"\nExample parent-child relationships:")
        for parent, children in list(extractor.parent_child_relations.items())[:5]:
            print(f"  {parent} -> {list(children)[:3]}...")
            
    except FileNotFoundError:
        print("Please download desc2020.xml from https://nlmpubs.nlm.nih.gov/projects/mesh/")
        print("and place it in the current directory")

if __name__ == "__main__":
    main()