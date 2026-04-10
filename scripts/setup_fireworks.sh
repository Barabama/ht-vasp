#!/bin/bash
# setup_fireworks.sh - Quick setup script for FireWorks + MongoDB integration
# Run this script on both nodes after updating configuration values

set -e

echo "=========================================="
echo "FireWorks + MongoDB Setup Script"
echo "=========================================="
echo ""

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

# Configuration (UPDATE THESE VALUES!)
MONGODB_HOST="192.168.9.144"  # Master node IP (429pro)
MONGODB_PORT="27017"
FIREWORKS_DB="fireworks_db"
FIREWORKS_USER="fireworks_user"
FIREWORKS_PASSWORD="fw123456"  # CHANGE THIS!
ATOMATE_DB="atomate_db"
ATOMATE_USER="atomate_user"
ATOMATE_PASSWORD="am123456"      # CHANGE THIS!

# Auto-detect node role based on hostname
# NOTE: Update this section if you have different hostnames or more nodes
CURRENT_HOSTNAME=$(hostname)
if [ "$CURRENT_HOSTNAME" = "429pro" ]; then
    NODE_ROLE="master"  # MongoDB + LaunchPad + Worker
elif [ "$CURRENT_HOSTNAME" = "429e" ]; then
    NODE_ROLE="compute"  # Worker only
else
    echo -e "${RED}Error: Unknown hostname '$CURRENT_HOSTNAME'${NC}"
    echo "Please update the script to recognize your hostname."
    echo ""
    echo "Example: Add your hostname to the detection logic:"
    echo "  elif [ \"\$CURRENT_HOSTNAME\" = \"your_hostname\" ]; then"
    echo "      NODE_ROLE=\"compute\""
    exit 1
fi

CONFIG_DIR="/home/mcmf429/nfs_hdd/2025/gaominliang/ht-vasp/configs"
LOG_DIR="/home/mcmf429/nfs_hdd/2025/gaominliang/ht-vasp/logs/fireworks"

# Check if running as root for certain operations
check_root() {
    if [ "$EUID" -ne 0 ]; then
        echo -e "${YELLOW}Warning: Some operations may require sudo privileges${NC}"
    fi
}

# Step 1: Install dependencies
install_dependencies() {
    echo ""
    echo "Step 1: Installing dependencies..."

    # Check and install system dependencies
    echo "Checking system dependencies..."
    MISSING_DEPS=()
    
    if ! command -v curl &> /dev/null; then
        MISSING_DEPS+=("curl")
    fi
    
    if ! command -v gpg &> /dev/null; then
        MISSING_DEPS+=("gnupg")
    fi
    
    if [ ${#MISSING_DEPS[@]} -ne 0 ]; then
        echo -e "${YELLOW}Installing missing system packages: ${MISSING_DEPS[*]}${NC}"
        sudo apt update -qq
        sudo apt install -y -qq "${MISSING_DEPS[@]}" > /dev/null 2>&1
        echo -e "${GREEN}✓ System dependencies installed${NC}"
    else
        echo -e "${GREEN}✓ All system dependencies present${NC}"
    fi

    # Activate conda environment
    source /opt/miniconda3/etc/profile.d/conda.sh
    conda activate htvasp || {
        echo -e "${RED}Error: htvasp conda environment not found${NC}"
        echo "Please create it first: conda create -n htvasp python=3.10"
        exit 1
    }

    # Install FireWorks
    pip install fireworks --quiet
    echo -e "${GREEN}✓ FireWorks installed${NC}"
}

# Step 2: Setup MongoDB (only on master node)
setup_mongodb() {
    if [ "$NODE_ROLE" = "master" ]; then
        echo ""
        echo "Step 2: Setting up MongoDB on node1..."

        # Check if MongoDB is already installed
        if ! command -v mongod &> /dev/null; then
            echo "Installing MongoDB 7.0..."

            # Import GPG key (force overwrite if exists)
            curl -fsSL https://www.mongodb.org/static/pgp/server-7.0.asc | \
                sudo gpg --dearmor --yes --output /usr/share/keyrings/mongodb-server-7.0.gpg

            # Add repository
            echo "deb [ signed-by=/usr/share/keyrings/mongodb-server-7.0.gpg ] \
http://repo.mongodb.org/apt/ubuntu jammy/mongodb-org/7.0 multiverse" | \
                sudo tee /etc/apt/sources.list.d/mongodb-org-7.0.list

            # Install
            sudo apt update
            sudo apt install -y mongodb-org

            echo -e "${GREEN}✓ MongoDB installed${NC}"
        else
            echo -e "${GREEN}✓ MongoDB already installed${NC}"
        fi

        # Configure MongoDB
        echo "Configuring MongoDB..."
        sudo cp /etc/mongod.conf /etc/mongod.conf.backup

        # Update bindIp to include localhost and the specified host
        # Use a more robust approach: replace the entire bindIp line
        if grep -q "^  bindIp:" /etc/mongod.conf; then
            sudo sed -i "s/^  bindIp:.*/  bindIp: 127.0.0.1,${MONGODB_HOST}/" /etc/mongod.conf
        elif grep -q "bindIp:" /etc/mongod.conf; then
            sudo sed -i "s/bindIp:.*/bindIp: 127.0.0.1,${MONGODB_HOST}/" /etc/mongod.conf
        else
            echo -e "${RED}Warning: Could not find bindIp in mongod.conf${NC}"
            echo "Adding bindIp configuration..."
            sudo sed -i '/^net:/a\  bindIp: 127.0.0.1,'"${MONGODB_HOST}" /etc/mongod.conf
        fi

        # Check if authorization is already enabled
        if ! grep -q "^security:" /etc/mongod.conf; then
            # Add security section at the end of file (no indentation)
            echo "" | sudo tee -a /etc/mongod.conf > /dev/null
            echo "security:" | sudo tee -a /etc/mongod.conf > /dev/null
            echo "  authorization: enabled" | sudo tee -a /etc/mongod.conf > /dev/null
        fi
        
        # Verify configuration by attempting to start MongoDB
        echo "Verifying MongoDB configuration..."
        sudo systemctl restart mongod
        sleep 3
        
        if sudo systemctl is-active --quiet mongod; then
            echo -e "${GREEN}✓ MongoDB configuration verified and service started${NC}"
        else
            echo -e "${RED}✗ MongoDB failed to start with new configuration${NC}"
            echo "Restoring backup..."
            sudo cp /etc/mongod.conf.backup /etc/mongod.conf
            sudo systemctl restart mongod
            echo ""
            echo "Error logs:"
            sudo journalctl -u mongod -n 20 --no-pager
            echo ""
            echo "Please check /etc/mongod.conf manually"
            echo "Common issues:"
            echo "  1. YAML indentation (security: must have no leading spaces)"
            echo "  2. Invalid option names or values"
            exit 1
        fi

        # Enable MongoDB service (already started during verification)
        sudo systemctl enable mongod

        echo -e "${GREEN}✓ MongoDB configured and started${NC}"
        
        # Wait for MongoDB to be ready
        echo "Waiting for MongoDB to be ready..."
        MAX_RETRIES=30
        RETRY_COUNT=0
        while ! mongosh --quiet --eval "db.adminCommand('ping')" &> /dev/null; do
            RETRY_COUNT=$((RETRY_COUNT + 1))
            if [ $RETRY_COUNT -ge $MAX_RETRIES ]; then
                echo -e "${RED}✗ MongoDB failed to start within timeout${NC}"
                echo "Check logs: sudo journalctl -u mongod -n 50"
                exit 1
            fi
            sleep 2
        done
        echo -e "${GREEN}✓ MongoDB is ready${NC}"

        # Create users
        echo ""
        echo "Creating MongoDB users..."
        
        # Check if users already exist
        if mongosh --quiet -u "${FIREWORKS_USER}" -p "${FIREWORKS_PASSWORD}" --authenticationDatabase "${FIREWORKS_DB}" "${FIREWORKS_DB}" --eval "db.getName()" &>/dev/null; then
            echo -e "${YELLOW}MongoDB users already exist, skipping creation${NC}"
        else
            echo "Creating admin and application users..."
            mongosh <<MONGOEOF
use admin
db.createUser({
    user: 'admin',
    pwd: '${FIREWORKS_PASSWORD}',
    roles: [ { role: 'userAdminAnyDatabase', db: 'admin' } ]
})

use ${FIREWORKS_DB}
db.createUser({
    user: '${FIREWORKS_USER}',
    pwd: '${FIREWORKS_PASSWORD}',
    roles: [
        { role: 'readWrite', db: '${FIREWORKS_DB}' },
        { role: 'dbAdmin', db: '${FIREWORKS_DB}' }
    ]
})

use ${ATOMATE_DB}
db.createUser({
    user: '${ATOMATE_USER}',
    pwd: '${ATOMATE_PASSWORD}',
    roles: [
        { role: 'readWrite', db: '${ATOMATE_DB}' },
        { role: 'dbAdmin', db: '${ATOMATE_DB}' }
    ]
})
MONGOEOF
            
            if [ $? -eq 0 ]; then
                echo -e "${GREEN}✓ MongoDB users created successfully${NC}"
            else
                echo -e "${RED}✗ Failed to create some MongoDB users${NC}"
                echo "Users may already exist. You can verify with:"
                echo "  mongosh -u admin -p <password> --authenticationDatabase admin"
            fi
        fi
    else
        echo ""
        echo "Step 2: Skipping MongoDB setup (not node1)..."
    fi
}

# Step 3: Configure FireWorks
configure_fireworks() {
    echo ""
    echo "Step 3: Configuring FireWorks..."

    # Create config directory
    mkdir -p ~/.fireworks
    mkdir -p "$LOG_DIR"
    mkdir -p "$CONFIG_DIR/fireworks"

    # Copy and customize mcmf_launchpad.yaml
    cp "$CONFIG_DIR/fireworks/mcmf_launchpad.yaml.template" ~/.fireworks/mcmf_launchpad.yaml
    sed -i "s/host: .*/host: ${MONGODB_HOST}/" ~/.fireworks/mcmf_launchpad.yaml
    sed -i "s/port: .*/port: ${MONGODB_PORT}/" ~/.fireworks/mcmf_launchpad.yaml
    sed -i "s/name: .*/name: ${FIREWORKS_DB}/" ~/.fireworks/mcmf_launchpad.yaml
    sed -i "s/username: .*/username: ${FIREWORKS_USER}/" ~/.fireworks/mcmf_launchpad.yaml
    sed -i "s/password: .*/password: ${FIREWORKS_PASSWORD}/" ~/.fireworks/mcmf_launchpad.yaml

    echo -e "${GREEN}✓ mcmf_launchpad.yaml configured${NC}"
    
    # Create symlink for FireWorks compatibility (FireWorks looks for my_launchpad.yaml by default)
    ln -sf ~/.fireworks/mcmf_launchpad.yaml ~/.fireworks/my_launchpad.yaml

    # Copy and customize mcmf_fworker.yaml
    if [ "$NODE_ROLE" = "master" ]; then
        cp "$CONFIG_DIR/fireworks/mcmf_fworker_node1.yaml.template" ~/.fireworks/mcmf_fworker.yaml
    else
        cp "$CONFIG_DIR/fireworks/mcmf_fworker_node2.yaml.template" ~/.fireworks/mcmf_fworker.yaml
    fi

    echo -e "${GREEN}✓ mcmf_fworker.yaml configured${NC}"
    
    # Create symlink for FireWorks compatibility
    ln -sf ~/.fireworks/mcmf_fworker.yaml ~/.fireworks/my_fworker.yaml

    # Copy mcmf_qadapter.yaml
    cp "$CONFIG_DIR/fireworks/mcmf_qadapter.yaml.template" ~/.fireworks/mcmf_qadapter.yaml

    echo -e "${GREEN}✓ mcmf_qadapter.yaml configured${NC}"
    
    # Create symlink for FireWorks compatibility
    ln -sf ~/.fireworks/mcmf_qadapter.yaml ~/.fireworks/my_qadapter.yaml
}

# Step 4: Configure atomate2 JobStore
configure_atomate2() {
    echo ""
    echo "Step 4: Configuring atomate2 JobStore..."

    # Copy and customize jobflow.yaml
    cp "$CONFIG_DIR/jobflow.yaml.template" "$CONFIG_DIR/jobflow.yaml"
    sed -i "s/host: .*/host: ${MONGODB_HOST}/g" "$CONFIG_DIR/jobflow.yaml"
    sed -i "s/port: .*/port: ${MONGODB_PORT}/g" "$CONFIG_DIR/jobflow.yaml"
    sed -i "s/username: atomate_user/username: ${ATOMATE_USER}/g" "$CONFIG_DIR/jobflow.yaml"
    sed -i "s/password: .*/password: ${ATOMATE_PASSWORD}/g" "$CONFIG_DIR/jobflow.yaml"

    echo -e "${GREEN}✓ jobflow.yaml configured${NC}"

    # Copy atomate2.yaml
    cp "$CONFIG_DIR/atomate2.yaml.template" "$CONFIG_DIR/atomate2.yaml"

    echo -e "${GREEN}✓ atomate2.yaml configured${NC}"

    # Set environment variables
    echo "" >> ~/.bashrc
    echo "# FireWorks and atomate2 configuration" >> ~/.bashrc
    echo "export JOBFLOW_CONFIG_FILE=$CONFIG_DIR/jobflow.yaml" >> ~/.bashrc
    echo "export ATOMATE2_CONFIG_FILE=$CONFIG_DIR/atomate2.yaml" >> ~/.bashrc

    export JOBFLOW_CONFIG_FILE="$CONFIG_DIR/jobflow.yaml"
    export ATOMATE2_CONFIG_FILE="$CONFIG_DIR/atomate2.yaml"

    echo -e "${GREEN}✓ Environment variables set${NC}"
}

# Step 5: Setup systemd service (optional)
setup_systemd_service() {
    echo ""
    echo "Step 5: Setting up systemd service (optional)..."
    read -p "Do you want to setup systemd service? (y/N): " setup_service

    if [ "$setup_service" = "y" ] || [ "$setup_service" = "Y" ]; then
        if [ "$NODE_ROLE" = "master" ]; then
            sudo cp "$CONFIG_DIR/systemd/fireworks-worker.service.template" \
                /etc/systemd/system/fireworks-worker.service
            sudo systemctl daemon-reload
            sudo systemctl enable fireworks-worker
            echo -e "${GREEN}✓ systemd service enabled (run 'sudo systemctl start fireworks-worker' to start)${NC}"
        else
            sudo cp "$CONFIG_DIR/systemd/fireworks-worker-node2.service.template" \
                /etc/systemd/system/fireworks-worker-node2.service
            sudo systemctl daemon-reload
            sudo systemctl enable fireworks-worker-node2
            echo -e "${GREEN}✓ systemd service enabled (run 'sudo systemctl start fireworks-worker-node2' to start)${NC}"
        fi
    else
        echo "Skipping systemd service setup."
        echo "You can manually start the worker with:"
        echo "  rlaunch -c $CONFIG_DIR/fireworks rapidfire"
    fi
}

# Step 6: Test connection
test_connection() {
    echo ""
    echo "Step 6: Testing connection..."

    # Test LaunchPad connection
    if lpad get_wflows &> /dev/null; then
        echo -e "${GREEN}✓ LaunchPad connection successful${NC}"
    else
        echo -e "${RED}✗ LaunchPad connection failed${NC}"
        echo "Please check your MongoDB credentials and network connectivity"
        exit 1
    fi

    # Show current status
    echo ""
    echo "Current FireWorks status:"
    lpad get_wflows
}

# Main execution
main() {
    echo "=========================================="
    echo "Detected Hostname: $CURRENT_HOSTNAME"
    echo "Node Role: $NODE_ROLE"
    if [ "$NODE_ROLE" = "master" ]; then
        echo "  - MongoDB Server: YES"
        echo "  - FireWorks LaunchPad: YES"
        echo "  - Worker: YES"
    else
        echo "  - MongoDB Server: NO (connect to $MONGODB_HOST)"
        echo "  - FireWorks LaunchPad: NO (use remote)"
        echo "  - Worker: YES"
    fi
    echo "=========================================="
    echo ""
    
    check_root
    install_dependencies
    setup_mongodb
    configure_fireworks
    configure_atomate2
    setup_systemd_service
    test_connection

    echo ""
    echo "=========================================="
    echo -e "${GREEN}Setup complete!${NC}"
    echo "=========================================="
    echo ""
    echo "Next steps:"
    echo "1. Review configuration files in ~/.fireworks/"
    echo "2. Test with a simple workflow:"
    echo "   python -c \"from fireworks import LaunchPad; lp = LaunchPad.auto_load(); print('Connected!')\""
    echo "3. Submit a test job:"
    echo "   lpad add workflows/test_workflow.yaml"
    echo "4. Start the worker (if not using systemd):"
    echo "   rlaunch -c $CONFIG_DIR/fireworks rapidfire"
    echo ""
    echo "Documentation: docs/fireworks-mongodb-deployment-guide.md"
}

# Run main function
main
