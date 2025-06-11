# MCP Server Troubleshooting Guide for UORCA

This guide helps you diagnose and fix issues with the MCP (Model Context Protocol) servers used in UORCA's AI-powered RNA-seq analysis features.

## 🚨 Quick Fix - Most Common Issues

If you're experiencing AI features not working in UORCA Explorer, try these steps in order:

### 1. Server Recovery (30 seconds)
```bash
cd UORCA/main_workflow/reporting
python tests/recover_mcp_servers.py --kill-all --restart
```

### 2. Refresh Your Browser
- If using Streamlit, refresh the page completely (Ctrl+F5 or Cmd+Shift+R)
- The app will reinitialize the MCP servers

### 3. Check Your Environment
```bash
python tests/test_mcp_servers.py --verbose
```

## 📋 Common Symptoms and Solutions

### Symptom: "MCP server is not running" Error
**What you see:** Error message in Streamlit about MCP servers not being available

**Quick fixes:**
1. **Kill stuck processes:**
   ```bash
   python tests/recover_mcp_servers.py --kill-all
   ```

2. **Restart servers:**
   ```bash
   python tests/recover_mcp_servers.py --restart --results-dir /path/to/your/results
   ```

3. **Refresh browser page** (important for Streamlit)

### Symptom: AI Features Show "Not Available"
**What you see:** Sidebar shows "🤖 AI Landing Page: Not Available" or "API Key Required"

**Solutions:**
1. **Set API Key:**
   ```bash
   export OPENAI_API_KEY=your_api_key_here
   ```
   Or create a `.env` file in the project root:
   ```
   OPENAI_API_KEY=your_api_key_here
   ```

2. **Install required packages:**
   ```bash
   pip install pydantic-ai openai
   # or with uv:
   uv add pydantic-ai openai
   ```

### Symptom: Servers Start but Become Unresponsive
**What you see:** "⚠️ MCP Servers: Some Issues" in sidebar

**Solutions:**
1. **Use the built-in recovery button** in Streamlit sidebar
2. **Or manually recover:**
   ```bash
   python tests/recover_mcp_servers.py --check-only
   python tests/recover_mcp_servers.py --kill-all --restart
   ```

### Symptom: "Event Loop" or "Asyncio" Errors
**What you see:** Technical errors mentioning event loops or asyncio

**Solutions:**
1. **Restart your Python session** (close and reopen terminal/IDE)
2. **Kill all Python processes** and start fresh:
   ```bash
   pkill -f "streamlit\|python.*mcp"
   ```
3. **Restart Streamlit:**
   ```bash
   streamlit run uorca_explorer.py
   ```

## 🔧 Diagnostic Tools

### Full System Test
```bash
cd UORCA/main_workflow/reporting
python tests/test_mcp_servers.py --results-dir /path/to/your/data --verbose
```

### Check What's Running
```bash
python tests/recover_mcp_servers.py --check-only
```

### Force Clean All Processes
```bash
python tests/recover_mcp_servers.py --kill-all --force
```

## 🛠️ Manual Recovery Procedures

### Complete Reset (Nuclear Option)
When nothing else works:

1. **Kill everything:**
   ```bash
   pkill -f "mcp_"
   pkill -f "streamlit"
   python tests/recover_mcp_servers.py --kill-all --force
   ```

2. **Wait and restart:**
   ```bash
   sleep 5
   streamlit run uorca_explorer.py
   ```

3. **Refresh browser completely**

### Server-Specific Recovery
If only one server type is problematic:

1. **Check which servers are problematic:**
   ```bash
   python tests/test_mcp_servers.py --verbose | grep "❌"
   ```

2. **Look for specific error patterns in logs**

3. **Restart individual components as needed**

### Environment Reset
If you suspect environment issues:

1. **Check Python and packages:**
   ```bash
   python --version  # Should be 3.9+
   pip list | grep -E "(pydantic-ai|openai|asyncio)"
   ```

2. **Recreate virtual environment if using one:**
   ```bash
   # With uv:
   uv sync --reinstall
   
   # With pip:
   pip install --force-reinstall -r requirements.txt
   ```

## 📊 Understanding Server Status

### Sidebar Status Indicators
- **🔗 MCP Servers: X Healthy** - All servers working
- **⚠️ MCP Servers: Some Issues** - Partial failure, recovery possible
- **🔗 MCP Servers: None Active** - Complete failure, restart needed

### Log Messages to Watch For
- `"Successfully set up X MCP servers"` - Good
- `"Server health check failed"` - Individual server problem
- `"MCP server process terminated"` - Process crashed
- `"Event loop already running"` - Asyncio conflict

## 🔍 Advanced Troubleshooting

### Server Process Investigation
```bash
# Find MCP processes
ps aux | grep mcp_

# Check if processes are responsive
python -c "
import psutil
for p in psutil.process_iter():
    if 'mcp_' in ' '.join(p.cmdline()):
        print(f'PID {p.pid}: {p.status()}')"
```

### Log Analysis
Enable verbose logging to see detailed information:
```bash
python tests/test_mcp_servers.py --verbose --results-dir /your/path 2>&1 | tee mcp_debug.log
```

Look for patterns in the log:
- Connection timeouts
- Permission errors
- Import failures
- Process startup failures

### Network/Port Issues
MCP servers use stdio communication, but if you suspect network issues:
```bash
# Check for unusual network activity
netstat -tlnp | grep python

# Check system resources
top -p $(pgrep -f mcp_)
```

## 🚫 Prevention Strategies

### Regular Maintenance
1. **Periodically clean up processes:**
   ```bash
   # Add to cron or run weekly
   python tests/recover_mcp_servers.py --kill-all
   ```

2. **Monitor resource usage** if running long analysis sessions

3. **Keep packages updated:**
   ```bash
   uv sync --upgrade  # or pip install --upgrade
   ```

### Best Practices
1. **Always refresh browser** after server errors
2. **Don't run multiple Streamlit instances** simultaneously
3. **Close unused browser tabs** running UORCA
4. **Set API keys in environment** rather than hardcoding
5. **Use the built-in recovery buttons** before manual intervention

### Environment Setup
1. **Use a clean Python environment:**
   ```bash
   # With uv (recommended):
   uv venv
   uv sync
   
   # Or with conda:
   conda create -n uorca python=3.11
   conda activate uorca
   ```

2. **Set environment variables properly:**
   ```bash
   # In .bashrc or .zshrc:
   export OPENAI_API_KEY=your_key
   export UORCA_DEFAULT_RESULTS_DIR=/path/to/your/results
   ```

## 🆘 When to Seek Help

Contact support or check GitHub issues if:

1. **Recovery tools consistently fail** after multiple attempts
2. **Environment tests pass but servers still won't start**
3. **You see unfamiliar error messages** not covered here
4. **Performance is severely degraded** even when servers start

### Information to Provide
When reporting issues, include:
1. Output of `python tests/test_mcp_servers.py --verbose`
2. Python version and OS
3. Full error messages from Streamlit
4. Whether the issue is reproducible

## 📚 Technical Details

### Architecture Overview
- **AutoStartManager**: Handles server lifecycle in Streamlit
- **ReportingAgent**: Manages AI analysis workflow
- **MCP Servers**: Data extraction and analysis services
- **Health Monitoring**: Automatic detection and recovery

### Server Communication
- Servers use **stdio transport** (stdin/stdout communication)
- **Asyncio** manages concurrent operations
- **Process management** handles server lifecycle
- **Health checks** verify server responsiveness

### Recovery Mechanisms
1. **Process-level monitoring** - Detects dead processes
2. **Communication timeouts** - Identifies unresponsive servers
3. **Automatic restart** - Attempts recovery without user intervention
4. **Manual recovery tools** - For complex failure scenarios

---

## 📞 Quick Reference Commands

```bash
# Emergency reset
python tests/recover_mcp_servers.py --kill-all --restart

# Status check
python tests/recover_mcp_servers.py --check-only

# Full diagnostic
python tests/test_mcp_servers.py --verbose

# Environment check
python tests/test_mcp_servers.py

# Force cleanup
python tests/recover_mcp_servers.py --kill-all --force
```

Remember: **When in doubt, refresh your browser!** Most Streamlit-related issues resolve with a page refresh after server recovery.