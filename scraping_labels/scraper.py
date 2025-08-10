import xml.etree.ElementTree as ET
import requests
# import gzip
from collections import defaultdict
import json

class MeSHHierarchyExtractor:
    def __init__(self):
        self.mesh_hierarchy = defaultdict(list)
        self.descriptor_to_tree = defaultdict(list)
        self.tree_to_descriptor = {}
        self.parent_child_relations = defaultdict(set)
        
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
        Parse MeSH XML and extract hierarchical relationships
        """
        if isinstance(xml_content, bytes):
            xml_content = xml_content.decode('utf-8')
            
        root = ET.fromstring(xml_content)
        
        for descriptor_record in root.findall('.//DescriptorRecord'):
            # Get descriptor UI and name
            descriptor_ui = descriptor_record.find('DescriptorUI').text
            descriptor_name = descriptor_record.find('DescriptorName/String').text
            
            # Get tree numbers (hierarchy positions)
            tree_number_list = descriptor_record.find('TreeNumberList')
            if tree_number_list is not None:
                for tree_number_elem in tree_number_list.findall('TreeNumber'):
                    tree_number = tree_number_elem.text
                    
                    # Store mappings
                    self.descriptor_to_tree[descriptor_name].append(tree_number)
                    self.tree_to_descriptor[tree_number] = descriptor_name
                    
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
        
        return similar_pairs
    
    def save_hierarchy_data(self, output_file="mesh_hierarchy.json"):
        """
        Save extracted hierarchy data to JSON file
        """
        hierarchy_data = {
            "descriptor_to_tree": dict(self.descriptor_to_tree),
            "tree_to_descriptor": self.tree_to_descriptor,
            "parent_child_relations": {k: list(v) for k, v in self.parent_child_relations.items()},
            "similar_pairs": self.get_label_hierarchy_for_micol()
        }
        
        with open(output_file, 'w') as f:
            json.dump(hierarchy_data, f, indent=2)
        
        print(f"Hierarchy data saved to {output_file}")
        return hierarchy_data

def main():
    # Example usage
    extractor = MeSHHierarchyExtractor()
    
    # Option 1: Download automatically (may not work due to file sizes)
    # xml_content = extractor.download_mesh_xml(2025)
    # save above xml to a file
    # if xml_content:
    #     with open('desc2025.xml', 'wb') as f:
    #         f.write(xml_content)
    
    # Option 2: Load from local file (recommended)
    # Download desc2020.xml manually and load it
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
        print(f"Total similar pairs for training: {len(hierarchy_data['similar_pairs'])}")
        
        # Show example relationships
        print(f"\nExample parent-child relationships:")
        for parent, children in list(extractor.parent_child_relations.items())[:5]:
            print(f"  {parent} -> {list(children)[:3]}...")
            
    except FileNotFoundError:
        print("Please download desc2020.xml from https://nlmpubs.nlm.nih.gov/projects/mesh/")
        print("and place it in the current directory")

if __name__ == "__main__":
    main()